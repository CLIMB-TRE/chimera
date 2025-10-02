#!/usr/bin/env python3

import os
import csv
import argparse
from pathlib import Path
import sys
import pybktree
import Levenshtein
from tqdm import tqdm


def taxon_string_to_taxon_name(
    name_tree: pybktree.BKTree, taxon_string: str, args: argparse.Namespace
) -> str | None:
    lineage_parts = taxon_string.split(";")

    lineage_names = [x[3:] for x in lineage_parts]

    for part in reversed(lineage_names):
        max_dist = len(part) // 8
        taxon_list = name_tree.find(part, max_dist)

        if taxon_list:
            if len(taxon_list) > 1:
                if args.verbose:
                    print(
                        f"Multiple names found for taxon {part}: {taxon_list}. Using first match - {sorted(taxon_list)[0]}",
                        file=sys.stderr,
                    )
            human_readable = sorted(taxon_list)[0][1]
            return human_readable

    return None


def run(args: argparse.Namespace):
    names_dmp = os.path.join(args.taxonomy, "names.dmp")

    tax_id_name_lookup = {}
    name_tax_id_lookup = {}
    print(f"Loading NCBI taxonomy from {names_dmp}", file=sys.stderr)
    with open(names_dmp, "r") as f:
        for line in f:
            fields = [i.lstrip() for i in line.split("\t|")]
            taxon_id, name, name_type = fields[0], fields[1], fields[3]
            if "scientific name" in name_type:
                tax_id_name_lookup[taxon_id] = name
                name_tax_id_lookup[name] = taxon_id

    print(f"Loaded {len(tax_id_name_lookup)} taxIDs from names.dmp", file=sys.stderr)

    name_tree = pybktree.BKTree(Levenshtein.distance, tax_id_name_lookup.values())

    print("Built BK-tree for taxon name lookup", file=sys.stderr)

    out_rows = []

    with open(args.sylph_taxonomy, "r") as infile:
        reader = csv.DictReader(
            infile,
            delimiter="\t",
            fieldnames=["contig_fname", "taxon_string"],
        )

        for row in tqdm(reader, desc="Processing rows"):
            taxon_name = taxon_string_to_taxon_name(
                name_tree, row["taxon_string"], args
            )
            if args.verbose:
                print(
                    f"Resolved taxon string {row['taxon_string']} to taxon name {taxon_name}",
                    file=sys.stderr,
                )
            if not taxon_name:
                if args.verbose:
                    print(
                        f"Warning: could not find taxon name for taxon string {row['taxon_string']} (contig file {row['contig_fname']})",
                        file=sys.stderr,
                    )
                print(
                    "\t".join(
                        (
                            row["contig_fname"],
                            row["taxon_string"],
                            "0",
                            "Unclassified",
                        )
                    ),
                    file=sys.stdout,
                )
                continue

            tax_id = name_tax_id_lookup.get(taxon_name)
            if tax_id:
                print(
                    "\t".join(
                        (
                            row["contig_fname"],
                            row["taxon_string"],
                            tax_id,
                            tax_id_name_lookup[tax_id],
                        )
                    ),
                    file=sys.stdout,
                )
            else:
                if args.verbose:
                    print(
                        f"ERROR: could not find taxID for taxon name {taxon_name} (contig file {row['contig_fname']})",
                        file=sys.stderr,
                    )
                    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Rationalise taxon strings from GTDB style to NCBI taxIDs"
    )
    parser.add_argument(
        "taxonomy",
        type=Path,
        help="Directory containing the NCBI taxonomy files (nodes.dmp, names.dmp)",
    )
    parser.add_argument(
        "sylph_taxonomy",
        type=Path,
        help="Path to the sylph taxonomy tsv file",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
