#!/usr/bin/env python3

import os
import csv
import argparse
from pathlib import Path
import sys
import pybktree
import Levenshtein
from tqdm import tqdm
from collections import namedtuple
import re

TAXON = namedtuple("TAXON", ["tax_id", "human_readable", "rank"])

RANK_LOOKUP = {
    "d": set(["superkingdom"]),
    "p": set(["phylum"]),
    "c": set(["class"]),
    "o": set(["order"]),
    "f": set(["family"]),
    "g": set(["genus"]),
    "s": set(["species", "strain", "subspecies", "species group", "species complex"]),
}


def taxon_levenshtein_distance(t1: TAXON, t2: TAXON) -> int:
    return Levenshtein.distance(t1.human_readable, t2.human_readable)


def taxon_string_to_taxon_name(
    name_tree: pybktree.BKTree,
    taxon_string: str,
    tally: dict,
    args: argparse.Namespace,
) -> TAXON | None:

    lineage_parts = taxon_string.split(";")

    for name, rank in reversed([(x[3:], x[0]) for x in lineage_parts]):
        if args.verbose:
            print(f"Looking up taxon name {name} with rank {rank}", file=sys.stderr)

        max_dist = len(name) // 8
        taxon_list = sorted(
            name_tree.find(TAXON(tax_id=None, human_readable=name, rank=rank), max_dist)
        )

        if taxon_list:
            if len(taxon_list) > 1:
                if args.verbose:
                    print(
                        f"Multiple names found for taxon {name}: {taxon_list}.",
                        file=sys.stderr,
                    )

            for taxon in taxon_list:
                if taxon[1].rank in RANK_LOOKUP[rank]:
                    if args.verbose:
                        print(
                            f"Selected taxon {taxon} based on rank {rank}",
                            file=sys.stderr,
                        )
                    tally["fuzzy matched"] += 1
                    return taxon[1]

    return None


def run(args: argparse.Namespace):
    names_dmp = os.path.join(args.taxonomy, "names.dmp")
    nodes_dmp = os.path.join(args.taxonomy, "nodes.dmp")

    parsed_taxonomy = {}
    tax_id_name_lookup = {}
    name_tax_id_lookup = {}
    rank_lookup = {}

    print(f"Loading NCBI taxonomy names from {names_dmp}", file=sys.stderr)
    with open(names_dmp, "r") as f:
        for line in f:
            fields = [i.lstrip() for i in line.split("\t|")]
            taxon_id, name, name_type = fields[0], fields[1], fields[3]
            if "scientific name" in name_type:
                tax_id_name_lookup[taxon_id] = name

    print(f"Loading NCBI taxonomy nodes from {nodes_dmp}", file=sys.stderr)
    with open(nodes_dmp, "r") as f:
        for line in f:
            fields = line.split("\t|\t")
            taxon_id, rank = fields[0], fields[2]
            rank_lookup[taxon_id] = rank

    for tax_id, name in tax_id_name_lookup.items():
        rank = rank_lookup.get(tax_id, "no_rank")
        parsed_taxonomy[tax_id] = TAXON(tax_id=tax_id, human_readable=name, rank=rank)

    print(f"Parsed {len(parsed_taxonomy)} taxa from NCBI taxonomy", file=sys.stderr)

    name_tree = pybktree.BKTree(taxon_levenshtein_distance, parsed_taxonomy.values())

    print("Built BK-tree for taxon name lookup", file=sys.stderr)

    out_rows = []

    with open(args.sylph_taxonomy, "r") as infile:
        reader = csv.DictReader(
            infile,
            delimiter="\t",
            fieldnames=["contig_fname", "taxon_string"],
        )

        sylph_taxonomy = []
        for row in tqdm(reader, desc="Pre-processing taxonomy strings"):
            subbed_taxon_string = re.sub(r"_[A-Za-z]\s", " ", row["taxon_string"])
            sylph_taxonomy.append(
                {
                    "contig_fname": row["contig_fname"],
                    "taxon_string": row["taxon_string"],
                    "subbed_taxon_string": subbed_taxon_string,
                }
            )

        assigned_taxa = {}
        tally = {"fuzzy matched": 0, "cached": 0, "not matched": 0}

        for i, row in enumerate(
            tqdm(sylph_taxonomy, desc="Matching taxon strings to NCBI taxIDs")
        ):
            if i > 0 and i % 10000 == 0:
                print(f"Processed {i} rows, current tally:\n{tally}", file=sys.stderr)

            taxon = assigned_taxa.get(row["subbed_taxon_string"])
            if not taxon:
                taxon = taxon_string_to_taxon_name(
                    name_tree,
                    row["subbed_taxon_string"],
                    tally,
                    args,
                )
                if taxon:
                    assigned_taxa[row["subbed_taxon_string"]] = taxon
            else:
                if args.verbose:
                    print(
                        f"Re-used cached taxon for taxon string {row['subbed_taxon_string']} -> {taxon}",
                        file=sys.stderr,
                    )
                tally["cached"] += 1

            if args.verbose:
                print(
                    f"Resolved taxon string {row['subbed_taxon_string']} to taxon name {taxon.human_readable if taxon else 'None'}",
                    file=sys.stderr,
                )
            if not taxon:
                if args.verbose:
                    print(
                        f"Warning: could not find taxon name for taxon string {row['taxon_string']} (contig file {row['contig_fname']})",
                        file=sys.stderr,
                    )
                tally["not matched"] += 1
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

            print(
                "\t".join(
                    (
                        row["contig_fname"],
                        row["taxon_string"],
                        taxon.tax_id,
                        taxon.human_readable,
                    )
                ),
                file=sys.stdout,
            )


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
