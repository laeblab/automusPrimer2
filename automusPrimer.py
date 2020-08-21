#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
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

RE_NUM_AMPLICONS = re.compile(rb"Descriptions of \[\s*(\d+)\s*\] potential amplicons")
RE_NUM_DIMERS = re.compile(rb"Dimer List \(\s*(\d+)\s*\)")
RE_NUM_HAIRPINS = re.compile(rb"Hairpin List \(\s*(\d+)\s*\)")


class TableError(RuntimeError):
    pass


def reverse_complement(seq):
    # Very simple (and inefficient) implementation
    return "".join(COMPLEMENT[nuc] for nuc in seq[::-1])


def safe_guide_rna_name(name):
    return "".join(char for char in name if char.isalnum() or char in " ._")


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


def run_mfeprimer3(params):
    log = logging.getLogger("run_mfeprimer3")
    args, candidate = params

    forward_primer = candidate["fwd_primer"]
    reverse_primer = candidate["rev_primer"]

    with tempfile.NamedTemporaryFile("wt") as handle:
        handle.write(f">forward\n{forward_primer}\n")
        handle.write(f">reverse\n{reverse_primer}\n")
        handle.flush()

        amplicons = run_and_parse_mfeprimer3(
            log,
            [args.mfeprimer3, "--db", args.fasta, "--in", handle.name],
            RE_NUM_AMPLICONS,
        )

        dimers = run_and_parse_mfeprimer3(
            log, [args.mfeprimer3, "dimer", "--in", handle.name], RE_NUM_DIMERS,
        )

        hairpins = run_and_parse_mfeprimer3(
            log, [args.mfeprimer3, "hairpin", "--in", handle.name], RE_NUM_HAIRPINS,
        )

    candidate["qc"] = (amplicons, dimers, hairpins)

    return candidate


def write_mfe3primer_cache(args, name, cache):
    filepath = args.output_folder / "mfeprimer3" / f"{name}.tsv"

    with filepath.open("wt") as handle:
        header = ["Forward Primer", "Reverse Primer", "Products", "Dimers", "Hairpins"]
        handle.write("\t".join(header) + "\n")

        for (fwd, rev), (prod, dimers, hairp) in sorted(cache.items()):
            handle.write("%s\t%s\t%s\t%s\t%s\n" % (fwd, rev, prod, dimers, hairp))


def read_mfe3primer_cache(args, name):
    cache = {}
    filepath = args.output_folder / "mfeprimer3" / f"{name}.tsv"
    filepath.parent.mkdir(parents=True, exist_ok=True)
    if filepath.exists():
        for row in read_table(filepath):
            forward_primer = row["Forward Primer"].upper()
            reverse_primer = row["Reverse Primer"].upper()

            cache[(forward_primer, reverse_primer)] = (
                int(row["Products"]),
                int(row["Dimers"]),
                int(row["Hairpins"]),
            )

    return cache


def find_best_primer_pairs(args, name, primer_pairs):
    log = logging.getLogger("find_best_primer_pairs")
    fwd_primers, rev_primers = primer_pairs

    n_fwd_primers = len(fwd_primers)
    n_rev_primers = len(rev_primers)

    max_depth = min(n_fwd_primers + n_rev_primers - 1, args.max_acceptable_depth)

    best_pair = None
    best_qc = (float("inf"), float("inf"), float("inf"))
    cache = read_mfe3primer_cache(args, name)

    candidates = []
    for depth in range(max_depth):
        for fwd_depth in range(depth + 1):
            if fwd_depth < n_fwd_primers and depth - fwd_depth < n_rev_primers:
                forward_primer = fwd_primers[fwd_depth]
                reverse_primer = rev_primers[depth - fwd_depth]

                cached_qc = cache.get((forward_primer, reverse_primer))
                if cached_qc is not None and (1, 0, 0) <= cached_qc <= best_qc:
                    best_pair = (forward_primer, reverse_primer, cached_qc)
                    best_qc = cached_qc

                if (forward_primer, reverse_primer) not in cache:
                    candidates.append(
                        (
                            args,
                            {
                                "fwd_primer": forward_primer,
                                "rev_primer": reverse_primer,
                                "qc": None,
                            },
                        )
                    )

    if best_qc == (1, 0, 0):
        log.info("Found cached f %s, r %s, qc %s", *best_pair)
        return best_pair
    elif best_pair is not None:
        log.info("Starting from f %s, r %s, qc %s", *best_pair)

    try:
        with multiprocessing.Pool() as pool:
            for candidate in pool.imap(run_mfeprimer3, candidates):
                if candidate is None:
                    break

                fwd_primer = candidate["fwd_primer"]
                rev_primer = candidate["rev_primer"]
                results = candidate["qc"]

                cache[(fwd_primer, rev_primer)] = results

                log.info("Tested f %s, r %s, qc %s", fwd_primer, rev_primer, results)
                if results == (1, 0, 0):
                    log.info("Found good primer-pair with single PCR product!")
                    return [fwd_primer, rev_primer, results]
                elif (1, 0, 0) <= results < best_qc:
                    best_pair = [fwd_primer, rev_primer, results]
                    best_qc = results

        if best_pair is None:
            log.error("No primer pairs found")
            return None

        log.warning("Found primer-pair QC %s", best_qc)
        return best_pair
    finally:
        write_mfe3primer_cache(args, name, cache)


def parse_primer3(args, name):
    def parse_primer3_file(handle):
        result = []
        for line in handle:
            fields = line.split()
            if fields and fields[0].isdigit():
                result.append(fields[1].upper())
        return result

    primer3_folder = args.output_folder / "primer3"

    with (primer3_folder / f"{name}.for").open("rt") as handle:
        forward_primers = parse_primer3_file(handle)

    with (primer3_folder / f"{name}.rev").open("rt") as handle:
        reverse_primers = parse_primer3_file(handle)

    return (forward_primers, reverse_primers)


def run_primer3(args, name, sequence):
    primer3_folder = args.output_folder / "primer3"
    primer3_folder.mkdir(parents=True, exist_ok=True)
    primer3_settings = primer3_folder / "settings"

    log = logging.getLogger("run_primer3")
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

    command = [args.primer3, "settings"]

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

    if not run_primer3(args, safe_name, target_region):
        return False

    suggested_primer_pairs = parse_primer3(args, safe_name)
    if suggested_primer_pairs is None:
        return False

    forward_primers, reverse_primers = suggested_primer_pairs
    forward_primers = [seq for seq in forward_primers if seq not in used_primers]
    reverse_primers = [seq for seq in reverse_primers if seq not in used_primers]
    suggested_primer_pairs = forward_primers, reverse_primers

    best_pair = find_best_primer_pairs(args, safe_name, suggested_primer_pairs)
    if best_pair is None:
        return False

    forward_primer, reverse_primer, qc = best_pair

    row["Forward Primer"] = forward_primer
    row["Reverse Primer"] = reverse_primer
    row["Quality Control"] = ":".join(str(value) for value in qc)

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

    # Write forward gRNA primers (Cas9)
    with (args.output_folder / "plate1_gRNAfwd.csv").open("wt") as handle:
        template = "%s,%s_gRNAfwd,ggaaaggacgaaacaccg%sgttttagagctagaaat\n"

        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(template % (pos, row["Safe Name"], row["Sequence"][1:-3]))

    # Write reverse gRNA primers (Cas9)
    with (args.output_folder / "plate2_gRNArev.csv").open("wt") as handle:
        template = "%s,%s_gRNArev,ctaaaac%scggtgtttcgtcctttccacaagatat\n"

        for pos, row in zip(PLATE_POSITIONS, table):
            seq = reverse_complement(row["Sequence"][1:-3])
            handle.write(template % (pos, row["Safe Name"], seq))

    # Write forward miSeq primers
    with (args.output_folder / "plate3_miseqfwd.csv").open("wt") as handle:
        template = "%s,%s_MiSeqfwd,tcgtcggcagcgtcagatgtgtataagagacag%s\n"

        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(template % (pos, row["Safe Name"], row["Forward Primer"]))

    # Write reverse miSeq primers
    with (args.output_folder / "plate4_miseqrev.csv").open("wt") as handle:
        template = "%s,%s_MiSeqrev,gtctcgtgggctcggagatgtgtataagagacag%s\n"

        for pos, row in zip(PLATE_POSITIONS, table):
            handle.write(template % (pos, row["Safe Name"], row["Reverse Primer"]))

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

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
