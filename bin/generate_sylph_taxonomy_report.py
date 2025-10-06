#!/usr/bin/env python3

import csv
from pathlib import Path
import sys


def genome_file_to_gcf_acc(file_name: str) -> str:
    if "ASM" in file_name:
        return file_name.split("/")[-1].split("_ASM")[0]
    return file_name.split("/")[-1].split("_genomic")[0]


def read_sylph_report(sylph_report_path: Path):
    """
    Reads a Sylph report file and returns a list of dictionaries representing each row.

    Args:
        sylph_report_path (str): Path to the Sylph report file.

    Returns:
        list: A list of dictionaries, each representing a row in the report.
    """

    with open(sylph_report_path, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        report_data = [row for row in reader]

    return report_data


def read_sylph_taxonomy(sylph_taxonomy_path: Path):
    """
    Reads a Sylph taxonomy file and returns a list of dictionaries representing each row.

    Args:
        sylph_taxonomy_path (str): Path to the modified Sylph taxonomy file. This file should have four columns: contig_fname, taxon_string, tax_id, and scientific name for the taxid.

    Returns:
        list: A list of dictionaries, each representing a row in the taxonomy.
    """

    with open(sylph_taxonomy_path, "r") as file:
        reader = csv.DictReader(
            file,
            delimiter="\t",
            fieldnames=["contig_fname", "taxon_string", "tax_id", "human_readable"],
        )
        taxonomy_data = {row["contig_fname"]: row for row in reader}

    return taxonomy_data


def add_taxon_data_to_report(report_data: list, taxonomy_data: dict) -> list:
    """
    Merges Sylph report data with taxonomy data based on contig filenames.

    Args:
        report_data (list): List of dictionaries representing the Sylph report.
        taxonomy_data (dict): Dictionary of dictionaries representing the Sylph taxonomy, keyed by contig_fname.

    Returns:
        list: A list of dictionaries, each representing a merged row from the report and taxonomy.
    """

    merged_data = []
    for row in report_data:
        contig_fname = genome_file_to_gcf_acc(row["contig_fname"])
        taxon_info = taxonomy_data.get(contig_fname)
        if not taxon_info:
            raise ValueError(
                f"Contig filename {contig_fname} not found in taxonomy data."
            )
        merged_row = {**row, **taxon_info}
        merged_data.append(merged_row)

    return merged_data


def run(args):
    report_data = read_sylph_report(args.sylph_report)
    taxonomy_data = read_sylph_taxonomy(args.sylph_taxonomy)
    merged_data = add_taxon_data_to_report(report_data, taxonomy_data)

    if not merged_data:
        sys.exit(2)

    # Output the merged data
    fieldnames = list(merged_data[0].keys())
    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in merged_data:
        writer.writerow(row)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a combined Sylph report with taxonomy information."
    )
    parser.add_argument(
        "--sylph_report",
        type=Path,
        required=True,
        help="Path to the Sylph report file.",
    )
    parser.add_argument(
        "--sylph_taxonomy",
        type=Path,
        required=True,
        help="Path to the modified Sylph taxonomy file.",
    )

    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
