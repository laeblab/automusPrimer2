#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import collections
import hashlib
import json
import logging
import multiprocessing
import re
import shutil
import subprocess
import sys
import tempfile

from pathlib import Path

import coloredlogs
import pysam

__version_tuple__ = (4, 0, 0)
__version__ = ".".join(map(str, __version_tuple__))


# Don't accept/cache amplicons if there are more than 3 for a pair of primers
MAX_AMPLICONS = 3

REQUIRED_HEADERS = (
    "Name",
    "Contig",
    "Sequence",
    "Position",
)

ADDED_HEADERS = ("Forward Primer", "Reverse Primer", "Quality Control", "PCR Sequence")

PLATE_POSITIONS = tuple(
    f"{row}{column}" for column in range(1, 17) for row in "ABCDEFGH"
)

COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "n": "n",
}

RE_NUM_DIMERS = re.compile(rb"Dimer List \(\s*(\d+)\s*\)")
RE_NUM_HAIRPINS = re.compile(rb"Hairpin List \(\s*(\d+)\s*\)")


class TableError(RuntimeError):
    pass


class Cas9GuideRNA:
    @classmethod
    def forward(cls, sequence):
        sequence = cls.process(sequence)

        return f"ggaaaggacgaaacacc{sequence}gttttagagctagaaat"

    @classmethod
    def reverse(cls, sequence):
        sequence = reverse_complement(cls.process(sequence))

        return f"ctaaaac{sequence}ggtgtttcgtcctttccacaagatat"

    @classmethod
    def process(cls, sequence):
        sequence = sequence.upper()
        assert len(sequence) == 23, repr(sequence)
        assert sequence.endswith("GG"), repr(sequence)

        sequence = "g" + sequence[1:-3]
        assert len(sequence) == 20, repr(sequence)

        return sequence


class MAD7GuideRNA:
    @classmethod
    def forward(cls, sequence):
        sequence = cls.process(sequence)

        return f"at{sequence}tttttttctagacccagct"

    @classmethod
    def reverse(cls, sequence):
        sequence = reverse_complement(cls.process(sequence))

        return f"agaaaaaaa{sequence}atctacaagagt"

    @classmethod
    def process(cls, sequence):
        sequence = sequence.upper()
        assert len(sequence) == 25, repr(sequence)
        assert sequence.startswith("TTT") or sequence.startswith("CTT"), repr(sequence)

        sequence = sequence[4:]
        assert len(sequence) == 21, repr(sequence)
        return sequence


class MiSeqPrimer:
    @classmethod
    def forward(cls, sequence):
        return f"tcgtcggcagcgtcagatgtgtataagagacag{sequence}"

    @classmethod
    def reverse(cls, sequence):
        return f"gtctcgtgggctcggagatgtgtataagagacag{sequence}"


def reverse_complement(seq):
    # Very simple (and inefficient) implementation
    return "".join(COMPLEMENT[nuc] for nuc in seq[::-1])


def safe_guide_rna_name(name):
    return "".join(char for char in name if char.isalnum() or char in " ._")


def hash_sequence(sequence):
    sequence = sequence.upper()

    bad_bases = set(sequence) - set("ACGT")
    if bad_bases:
        return None

    hasher = hashlib.new("md5")
    hasher.update(sequence.encode("utf8"))

    return hasher.hexdigest().upper()


def log_output(log, label, data):
    text = data.decode("utf-8", errors="replace")
    for line in text.split("\n"):
        line = line.rstrip()
        if line:
            log("%s%s", label, line)


def run_command(log, label, command, cwd=None):
    log.debug("Running command %s", [str(value) for value in command])
    process = subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=cwd,
    )

    stdout, stderr = process.communicate()
    log_output(log.debug, f"[{label}] ", stdout)
    log_output(log.error, f"[{label}] ", stderr)

    if process.returncode:
        log.error("Failed to run %s (returncode = %i)", label, process.returncode)
        return False

    return True


def read_table(filepath):
    with filepath.open("rt") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if not header:
            raise TableError(f"No header found in '{filepath}'")

        for linenum, line in enumerate(handle, start=2):
            if line.strip() and not line.startswith("#"):
                fields = [value.strip() for value in line.split("\t")]
                if len(fields) != len(header):
                    raise TableError(
                        f"Wrong number of columns on {linenum}, "
                        f"expected {len(header)} columns, but found {len(fields)}"
                    )

                yield dict(zip(header, fields))


def run_and_parse_mfeprimer3(log, command, regexp):
    log.debug("Running command %s", [str(value) for value in command])
    process = subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = process.communicate()
    log_output(log.error, "[mfeprimer3] ", stderr)

    if process.returncode:
        log.error("Failed to run mfeprimer3 (returncode = %i)", process.returncode)
        None

    match = regexp.search(stdout)
    if match is None:
        log.error("Regexp not found for command %r", command)
        return None

    return int(match.groups()[0])


def run_and_collect_mfeprimer3_amplicons(log, args, fastapath, outpath):
    command = [
        args.mfeprimer3,
        "--db",
        args.fasta,
        "--in",
        fastapath,
        "--json",
        "--out",
        outpath / "out",
    ]

    log.debug("Running command %s", [str(value) for value in command])
    process = subprocess.Popen(
        command,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = process.communicate()
    log_output(log.debug, "[mfeprimer3] ", stdout)
    log_output(log.error, "[mfeprimer3] ", stderr)

    if process.wait():
        log.error("Failed to run mfeprimer3 (returncode = %i)", process.returncode)
        return None

    with (outpath / "out.json").open() as handle:
        data = json.load(handle)

    # This seems to happen when mfeprimer3 hits memory limits
    if data["AmpList"] is None:
        log.error(
            "mfeprimer3 did not produce amplicon list (returncode = %i)",
            process.returncode,
        )
        return None

    # GuideRNAs may target duplicated regions / regions with alternative scaffolds;
    # in that case we do not wish to exclude primer-pairs that amplify both regions.
    # TODO: Handle very similar (identical length) amplicons
    # TODO: Prioritize globally unique amplicons if possible
    amplicons = set()
    for amplicon in data["AmpList"]:
        amplicons.add(amplicon["P"]["Seq"]["Seq"].upper())

    return tuple(amplicons)


def run_mfeprimer3(params):
    log = logging.getLogger("run_mfeprimer3")
    args, forward_primer, reverse_primer = params

    with tempfile.TemporaryDirectory() as tempdir:
        temppath = Path(tempdir)

        fastapath = temppath / "in.fasta"
        with fastapath.open("wt") as handle:
            handle.write(f">forward\n{forward_primer}\n")
            handle.write(f">reverse\n{reverse_primer}\n")

        amplicons = run_and_collect_mfeprimer3_amplicons(
            log=log, args=args, fastapath=fastapath, outpath=temppath
        )

        dimers = run_and_parse_mfeprimer3(
            log, [args.mfeprimer3, "dimer", "--in", fastapath], RE_NUM_DIMERS,
        )

        hairpins = run_and_parse_mfeprimer3(
            log, [args.mfeprimer3, "hairpin", "--in", fastapath], RE_NUM_HAIRPINS,
        )

    if None in (amplicons, dimers, hairpins):
        return None

    return {
        "fwd_primer": forward_primer,
        "rev_primer": reverse_primer,
        "qc": [len(amplicons), dimers, hairpins],
        "amplicons": amplicons,
    }


def write_mfe3primer_cache(args, name, cache):
    filepath = args.output_folder / "mfeprimer3" / f"{name}.json"

    with filepath.open("wt") as handle:
        json.dump(list(cache.items()), handle)


def read_mfe3primer_cache(args, name):
    cache = collections.OrderedDict()
    filepath = args.output_folder / "mfeprimer3" / f"{name}.json"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    if filepath.exists():
        with filepath.open("rt") as handle:
            for key, values in json.load(handle):
                cache[tuple(key)] = values

    return cache


def primer_permutations(args, primer_pairs):
    fwd_primers, rev_primers = primer_pairs

    n_fwd_primers = len(fwd_primers)
    n_rev_primers = len(rev_primers)

    max_depth = min(n_fwd_primers + n_rev_primers - 1, args.max_acceptable_depth)

    for depth in range(max_depth):
        for fwd_depth in range(depth + 1):
            if fwd_depth < n_fwd_primers and depth - fwd_depth < n_rev_primers:
                forward_primer = fwd_primers[fwd_depth]
                reverse_primer = rev_primers[depth - fwd_depth]

                yield forward_primer, reverse_primer


def find_best_primer_pairs(args, name, primer_pairs):
    log = logging.getLogger("find_best_primer_pairs")
    cache = read_mfe3primer_cache(args, name)
    candidates = []

    best_pair = None
    best_qc = [float("inf"), float("inf"), float("inf")]
    for forward_primer, reverse_primer in primer_permutations(args, primer_pairs):
        cached = cache.get((forward_primer, reverse_primer))
        if cached is not None:
            if [1, 0, 0] <= cached["qc"] < best_qc:
                best_qc = cached["qc"]
                best_pair = {
                    "forward": forward_primer,
                    "reverse": reverse_primer,
                    "qc": cached["qc"],
                    "amplicons": cached["amplicons"],
                }
        else:
            candidates.append((args, forward_primer, reverse_primer))

    if best_qc == [1, 0, 0]:
        log.info("Found cached f %(forward)s, r %(reverse)s, qc %(qc)s" % best_pair)
        return best_pair
    elif best_pair is not None:
        log.info("Starting from f %(forward)s, r %(reverse)s, qc %(qc)s" % best_pair)

    try:
        with multiprocessing.Pool(args.threads) as pool:
            for candidate in pool.imap(run_mfeprimer3, candidates):
                if candidate is None:
                    break

                fwd_primer = candidate["fwd_primer"]
                rev_primer = candidate["rev_primer"]
                results = candidate["qc"]
                amplicons = candidate["amplicons"]
                if len(amplicons) > MAX_AMPLICONS:
                    amplicons = "N" * len(amplicons)

                cache[(fwd_primer, rev_primer)] = {
                    "qc": results,
                    "amplicons": amplicons,
                }

                log.info("Tested f %s, r %s, qc %s", fwd_primer, rev_primer, results)
                if [1, 0, 0] <= results < best_qc:
                    best_pair = {
                        "forward": fwd_primer,
                        "reverse": rev_primer,
                        "qc": results,
                        "amplicons": amplicons,
                    }
                    best_qc = results

                if results == [1, 0, 0]:
                    log.info("Found good primer-pair with single PCR product!")
                    return best_pair

        if best_pair is None:
            log.error("No primer pairs found")
            return None
        elif best_qc > [MAX_AMPLICONS, 0, 0]:
            log.error("Best primer pair produces too many amplicons %r", best_qc)
            return None

        log.warning("Found primer-pair QC %s", best_qc)
        return best_pair
    finally:
        write_mfe3primer_cache(args, name, cache)


def parse_primer3(args, sequence):
    def parse_primer3_file(handle):
        result = []
        for line in handle:
            fields = line.split()
            if fields and fields[0].isdigit():
                result.append(fields[1].upper())
        return result

    name = hash_sequence(sequence)
    primer3_folder = args.output_folder / "primer3"

    with (primer3_folder / f"{name}.for").open("rt") as handle:
        forward_primers = parse_primer3_file(handle)

    with (primer3_folder / f"{name}.rev").open("rt") as handle:
        reverse_primers = parse_primer3_file(handle)

    return (forward_primers, reverse_primers)


def run_primer3(args, seqname, sequence):
    log = logging.getLogger("run_primer3")

    name = hash_sequence(sequence)
    if name is None:
        log.error("Invalid sequence for %r; contains non-ACGT bases!", seqname)
        return False

    primer3_folder = args.output_folder / "primer3"
    primer3_folder.mkdir(parents=True, exist_ok=True)
    primer3_settings = primer3_folder / f"{name}.settings"

    output_for = primer3_settings.with_suffix(".for")
    output_rev = primer3_settings.with_suffix(".rev")
    if output_for.exists() and output_rev.exists():
        log.info("Primer3 results for %s already exist (%s)", seqname, name)
        return True

    log.info("Writing primer3 settings to '%s'", primer3_settings)
    with primer3_settings.open("wt") as handle:
        settings = [
            f"SEQUENCE_ID={name}",
            f"SEQUENCE_TEMPLATE={sequence}",
            # left primer begin, length, right primer begin, length
            "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,99,200,99",
            "PRIMER_TASK=pick_primer_list",
            "PRIMER_PRODUCT_SIZE_RANGE=100-250",
            "P3_FILE_FLAG=1",
            "=",
        ]

        for line in settings:
            print(line, file=handle)

    command = [args.primer3, primer3_settings.name]

    return run_command(log, "primer3", command, cwd=primer3_folder)


def pick_primers_for_guide_rna(args, fasta, row, used_primers):
    name = row["Name"]
    contig = row["Contig"]
    position = int(row["Position"])

    log = logging.getLogger("pick_primers_for_guide_rna")
    log.info("Finding primer for gRNA %r at %s:%s", name, contig, position)

    try:
        contig_len = fasta.get_reference_length(contig)
    except KeyError:
        log.error("Contig %r not found in FASTA file", contig)
        return False

    region_start = position - args.max_primer_distance
    if region_start < 0:
        log.warning("gRNA is too close to contig 5p: %r", name)
        region_start = 0

    region_end = position + args.max_primer_distance
    if region_end > contig_len:
        log.warning("gRNA is too close to contig 3p: %r", name)
        region_end = contig_len

    safe_name = row["Safe Name"]
    target_region = fasta.fetch(contig, region_start, region_end)
    if "N" in target_region.upper():
        log.error("Uncalled bases in region around %r", name)
        return False

    if not run_primer3(args, safe_name, target_region):
        return False

    suggested_primer_pairs = parse_primer3(args, target_region)
    if suggested_primer_pairs is None:
        return False

    forward_primers, reverse_primers = suggested_primer_pairs
    forward_primers = [seq for seq in forward_primers if seq not in used_primers]
    reverse_primers = [seq for seq in reverse_primers if seq not in used_primers]
    suggested_primer_pairs = forward_primers, reverse_primers

    result = find_best_primer_pairs(args, safe_name, suggested_primer_pairs)
    if result is None:
        return False

    forward_primer = result["forward"]
    reverse_primer = result["reverse"]

    row["Forward Primer"] = forward_primer
    row["Reverse Primer"] = reverse_primer
    row["Quality Control"] = ":".join(str(value) for value in result["qc"])
    row["Predicted Amplicons"] = result["amplicons"]

    used_primers.add(forward_primer)
    used_primers.add(reverse_primer)

    pcr_product = target_region.upper()
    pcr_product = pcr_product[pcr_product.index(forward_primer) :]
    pcr_product = pcr_product[
        : pcr_product.index(reverse_complement(reverse_primer)) + len(reverse_primer)
    ]

    row["PCR Sequence"] = pcr_product

    return True


def read_guide_rna_table(filepath):
    log = logging.getLogger("read_guide_rna_table")
    log.info("Reading gRNA table from '%s'", filepath)

    try:
        table = []
        for rownum, row in enumerate(read_table(filepath)):
            missing = sorted(set(REQUIRED_HEADERS) - set(row))
            if missing:
                log.error("Required column(s) not found in '%s': %s", filepath, missing)
                return None

            missing = []
            for key in REQUIRED_HEADERS:
                value = row[key]
                if not value:
                    missing.append(key)

            if missing:
                log.error("Missing values for %s on row %i", missing, rownum)
                return None

            for key in ADDED_HEADERS:
                row.setdefault(key, "")

            row["Safe Name"] = safe_guide_rna_name(row["Name"])

            table.append(row)

        return table
    except TableError as error:
        log.error("%s", error)
        return None


def setup_logging(level=logging.INFO):
    coloredlogs.install(
        level=level,
        fmt="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        level_styles={
            "info": {"color": "white"},
            "critical": {"color": "red", "bold": True},
            "verbose": {"color": "blue"},
            "error": {"color": "red"},
            "debug": {"color": "green"},
            "warning": {"color": "yellow"},
        },
    )


def parse_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("fasta", type=Path, help="MFEPrimer3 indexed FASTA file")
    parser.add_argument("table", type=Path)
    parser.add_argument("output_folder", type=Path)
    parser.add_argument(
        "--threads",
        type=lambda value: max(1, int(value)),
        default=multiprocessing.cpu_count(),
        help="Number of worker threads used while evaluating candidate primer pairs",
    )
    parser.add_argument(
        "--max-acceptable-depth",
        type=int,
        default=25,
        help="primer3 finds a number of fwd and rev primers. These are ranked. "
        "The depth is the sum of how far down these ranked lists we had to go",
    )
    parser.add_argument(
        "--max-primer-distance",
        type=int,
        default=150,
        help="Primer must be located within this many bp of cut-site",
    )
    parser.add_argument(
        "--log-level",
        type=str.upper,
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        default="INFO",
        help="Log level for output written to STDERR",
    )
    parser.add_argument(
        "--enzyme",
        type=str.upper,
        choices=("CAS9", "MAD7"),
        default="CAS9",
        help="Create gRNA constructs for the specified enzyme",
    )

    group = parser.add_argument_group("Executables")
    group.add_argument(
        "--primer3",
        type=Path,
        default=Path("primer3_core"),
        help="Location of 'primer3' executable",
    )
    parser.add_argument(
        "--mfeprimer3",
        type=Path,
        default=Path("mfeprimer3"),
        help="Location of 'mfeprimer3' executable",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    setup_logging(args.log_level)

    log = logging.getLogger("main")
    if not args.fasta.is_file():
        log.error("FASTA file not found at '%s'", args.fasta)
        return 1
    log.info("FASTA file found at '%s'", args.fasta)

    mfeprimer_db = args.fasta.with_suffix(args.fasta.suffix + ".primerqc")
    if not mfeprimer_db.is_file():
        log.error("MFEPrimer3 index not found at '%s'", mfeprimer_db)
        log.error('Please run \'mfeprimer3 index -i "%s"', args.fasta)
        return 1
    log.info("MFEPrimer3 index found at '%s'", mfeprimer_db)

    primer3 = shutil.which(args.primer3)
    if not primer3:
        log.error("primer3 executable not found at '%s'", args.primer3)
        return 1
    args.primer3 = primer3

    mfeprimer3 = shutil.which(args.mfeprimer3)
    if not mfeprimer3:
        log.error("mfeprimer3 executable not found at '%s'", args.mfeprimer3)
        return 1
    args.mfeprimer3 = mfeprimer3

    args.output_folder.mkdir(parents=True, exist_ok=True)

    table = read_guide_rna_table(args.table)
    if table is None:
        return 1

    if len(table) > len(PLATE_POSITIONS):
        log.error(
            "Table can at most contain %i rows, but contains %i rows",
            len(PLATE_POSITIONS),
            len(table),
        )
        return 1

    with pysam.FastaFile(args.fasta) as fasta:
        used_primers = set()
        for row in table:
            if not pick_primers_for_guide_rna(args, fasta, row, used_primers):
                return 1

    if args.enzyme == "CAS9":
        enzyme = Cas9GuideRNA
    elif args.enzyme == "MAD7":
        enzyme = MAD7GuideRNA
    else:
        raise NotImplementedError(f"Unknown enzyme {args.enzyme!r}")

    # Write forward gRNA primers (Cas9)
    with (args.output_folder / "plate1_gRNAfwd.csv").open("wt") as handle:
        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(
                "{},{}_gRNAfwd,{}\n".format(
                    pos, row["Safe Name"], enzyme.forward(row["Sequence"])
                )
            )

    # Write reverse gRNA primers (Cas9)
    with (args.output_folder / "plate2_gRNArev.csv").open("wt") as handle:
        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(
                "{},{}_gRNArev,{}\n".format(
                    pos, row["Safe Name"], enzyme.reverse(row["Sequence"])
                )
            )

    # Write forward miSeq primers
    with (args.output_folder / "plate3_miseqfwd.csv").open("wt") as handle:
        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(
                "{},{}_MiSeqfwd,{}\n".format(
                    pos, row["Safe Name"], MiSeqPrimer.forward(row["Forward Primer"]),
                )
            )

    # Write reverse miSeq primers
    with (args.output_folder / "plate4_miseqrev.csv").open("wt") as handle:
        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(
                "{},{}_MiSeqrev,{}\n".format(
                    pos, row["Safe Name"], MiSeqPrimer.reverse(row["Reverse Primer"]),
                )
            )

    # Write PCR products
    with (args.output_folder / "targets.fasta").open("wt") as handle:
        for row in table:
            handle.write(">%s\n%s\n" % (row["Safe Name"], row["PCR Sequence"]))

    # Write StyrKO table
    with (args.output_folder / "styrko.tsv").open("wt") as handle:
        header = REQUIRED_HEADERS + ADDED_HEADERS
        handle.write("%s\n" % ("\t".join(header),))

        for row in table:
            row = dict(row)
            row["Reverse Primer"] = reverse_complement(row["Reverse Primer"])
            row = [row[key] for key in header]
            handle.write("%s\n" % ("\t".join(row),))

    # Write StyrKO table
    with (args.output_folder / "predicted_amplicons.tsv").open("wt") as handle:
        handle.write("Name\tNth\tSequence\n")

        for row in table:
            for idx, sequence in enumerate(sorted(row["Predicted Amplicons"])):
                handle.write("%s\t%s\t%s\n" % (row["Name"], idx + 1, sequence))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
