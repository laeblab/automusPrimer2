#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import sys

from pathlib import Path

import openpyxl


PLATE_POSITIONS = tuple(
    f"{row}{column}" for column in range(1, 17) for row in "ABCDEFGH"
)

PLATES = [
    "plate1_gRNAfwd.csv",
    "plate2_gRNArev.csv",
    "plate3_miseqfwd.csv",
    "plate4_miseqrev.csv",
]


def write_plate(filepath, rows):
    workbook = openpyxl.Workbook()
    worksheet = workbook.active

    for idx, row in enumerate(rows):
        worksheet.append((PLATE_POSITIONS[idx], row["Name"], row["Sequence"]))

    workbook.save(str(filepath))


def read_plate(filepath):
    rows = []
    header = ["Position", "Name", "Sequence"]
    with filepath.open("rt") as handle:
        for line in handle:
            line = line.strip()
            if line:
                fields = line.split(",")
                if len(fields) != len(header):
                    raise RuntimeError(f"Wrong number of columns in '{filepath}'")

                rows.append(dict(zip(header, fields)))

    return rows


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("folders", nargs="+", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)

    for plate in PLATES:
        rows = []
        for folder in args.folders:
            rows.extend(read_plate(folder / plate))

        write_plate((args.output / plate).with_suffix(".xlsx"), rows)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
