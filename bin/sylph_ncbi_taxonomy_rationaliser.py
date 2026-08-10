#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import re
import sys
from collections import namedtuple
from pathlib import Path

TAXON = namedtuple("TAXON", ["tax_id", "human_readable", "rank"])

GTDB_ACCESSION_RE = re.compile(r"^(?:RS_|GB_)?(GC[AF]_\d+)(?:\.\d+)?$")

ALIAS_NAME_CLASSES = {
    "scientific name",
    "synonym",
    "genbank synonym",
    "equivalent name",
    "common name",
    "genbank common name",
    "includes",
}


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def normalise_accession(accession: str) -> tuple[str, str]:
    """Return (versioned, unversioned) normalised GCA/GCF accession, stripping RS_/GB_ prefixes."""
    match = GTDB_ACCESSION_RE.match(accession.strip())
    if not match:
        stripped = re.sub(r"^(RS_|GB_)", "", accession.strip())
        base = stripped.split(".")[0]
        return stripped, base
    versioned = accession.strip()
    versioned = re.sub(r"^(RS_|GB_)", "", versioned)
    return versioned, match.group(1)


def load_gtdb_metadata(paths: list[Path]) -> dict:
    """
    Load one or more GTDB metadata TSVs (optionally gzipped) into a lookup dict
    keyed by normalised accession (both versioned and version-stripped) -> row dict
    with the columns we need.
    """
    gtdb_lookup = {}
    wanted_cols = [
        "accession",
        "ncbi_genbank_assembly_accession",
        "ncbi_taxid",
        "ncbi_species_taxid",
        "ncbi_organism_name",
        "gtdb_taxonomy",
        "ncbi_taxonomy",
    ]

    for path in paths:
        print(f"Loading GTDB metadata from {path}", file=sys.stderr)
        with open_maybe_gzip(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                record = {col: row.get(col, "") for col in wanted_cols}
                for acc_col in ("accession", "ncbi_genbank_assembly_accession"):
                    acc = record.get(acc_col)
                    if not acc:
                        continue
                    versioned, unversioned = normalise_accession(acc)
                    gtdb_lookup.setdefault(versioned, record)
                    gtdb_lookup.setdefault(unversioned, record)

    print(f"Loaded {len(gtdb_lookup)} GTDB accession keys", file=sys.stderr)
    return gtdb_lookup


def lookup_gtdb_record(gtdb_lookup: dict, contig_fname: str) -> dict | None:
    versioned, unversioned = normalise_accession(contig_fname)
    return gtdb_lookup.get(versioned) or gtdb_lookup.get(unversioned)


def load_ncbi_taxonomy(taxonomy_dir: Path):
    names_dmp = os.path.join(taxonomy_dir, "names.dmp")
    nodes_dmp = os.path.join(taxonomy_dir, "nodes.dmp")
    merged_dmp = os.path.join(taxonomy_dir, "merged.dmp")
    delnodes_dmp = os.path.join(taxonomy_dir, "delnodes.dmp")

    parsed_taxonomy = {}
    tax_id_name_lookup = {}
    rank_lookup = {}
    name_alias_lookup = {}

    print(f"Loading NCBI taxonomy names from {names_dmp}", file=sys.stderr)
    with open(names_dmp, "r") as f:
        for line in f:
            fields = [i.lstrip() for i in line.split("\t|")]
            taxon_id, name, name_class = fields[0], fields[1], fields[3].strip()
            if name_class == "scientific name":
                tax_id_name_lookup[taxon_id] = name
            if name_class in ALIAS_NAME_CLASSES:
                name_alias_lookup.setdefault(name.lower(), set()).add(taxon_id)

    print(f"Loading NCBI taxonomy nodes from {nodes_dmp}", file=sys.stderr)
    with open(nodes_dmp, "r") as f:
        for line in f:
            fields = line.split("\t|\t")
            taxon_id, rank = fields[0], fields[2]
            rank_lookup[taxon_id] = rank

    for tax_id, name in tax_id_name_lookup.items():
        rank = rank_lookup.get(tax_id, "no_rank")
        parsed_taxonomy[tax_id] = TAXON(tax_id=tax_id, human_readable=name, rank=rank)

    merged_map = {}
    if os.path.isfile(merged_dmp):
        print(f"Loading merged taxIDs from {merged_dmp}", file=sys.stderr)
        with open(merged_dmp, "r") as f:
            for line in f:
                fields = [i.strip() for i in line.split("\t|")]
                merged_map[fields[0]] = fields[1]

    deleted_ids = set()
    if os.path.isfile(delnodes_dmp):
        print(f"Loading deleted taxIDs from {delnodes_dmp}", file=sys.stderr)
        with open(delnodes_dmp, "r") as f:
            for line in f:
                deleted_ids.add(line.split("\t|")[0].strip())

    print(f"Parsed {len(parsed_taxonomy)} taxa from NCBI taxonomy", file=sys.stderr)

    return parsed_taxonomy, name_alias_lookup, merged_map, deleted_ids


def resolve_tax_id(tax_id: str, parsed_taxonomy: dict, merged_map: dict, deleted_ids: set):
    """
    Check whether tax_id is present in the provided NCBI taxonomy dump.
    Returns (resolved_tax_id_or_None, status) where status is one of:
    "valid", "merged", "deleted", "absent".
    """
    if not tax_id:
        return None, "absent"
    if tax_id in parsed_taxonomy:
        return tax_id, "valid"
    if tax_id in merged_map:
        new_id = merged_map[tax_id]
        if new_id in parsed_taxonomy:
            return new_id, "merged"
    if tax_id in deleted_ids:
        return None, "deleted"
    return None, "absent"


def literal_name_match(name: str, name_alias_lookup: dict) -> list[str]:
    if not name:
        return []
    return sorted(name_alias_lookup.get(name.strip().lower(), set()))


def strip_gtdb_suffix(name: str) -> str:
    return re.sub(r"_[A-Z]{1,2}(?=[\s;]|$)", "", name)


def species_name_from_taxon_string(taxon_string: str) -> str | None:
    for part in reversed(taxon_string.split(";")):
        if part.startswith("s__") and part[3:]:
            return strip_gtdb_suffix(part[3:])
    return None


class ReportWriter:
    FIELDNAMES = [
        "contig_fname",
        "taxon_string",
        "gtdb_accession",
        "gtdb_ncbi_taxid",
        "gtdb_ncbi_organism_name",
        "gtdb_taxonomy",
        "failure_reason",
        "suggested_tax_id",
        "suggested_name",
        "accepted",
        "resolution_method",
    ]

    def __init__(self, path):
        self._file = open(path, "w", newline="")
        self._writer = csv.DictWriter(
            self._file, fieldnames=self.FIELDNAMES, delimiter="\t", lineterminator="\n"
        )
        self._writer.writeheader()
        self._file.flush()

    def write(self, row: dict):
        self._writer.writerow(row)
        self._file.flush()

    def close(self):
        self._file.close()


class Tty:
    """Lazily opens /dev/tty so a run with no unresolved taxa needs no controlling terminal."""

    def __init__(self):
        self._handle = None

    def readline(self):
        if self._handle is None:
            self._handle = open("/dev/tty", "r")
        return self._handle.readline()

    def close(self):
        if self._handle is not None:
            self._handle.close()


def confirm_match(tty, taxon_string: str, contig_fname: str, gtdb_record: dict | None, suggestions: list) -> tuple[str, str] | None:
    """
    Prompt the user (via /dev/tty) to confirm a suggested tax_id/name match.
    Returns (tax_id, human_readable) if accepted, else None.
    """
    print("\n--- Unresolved taxon ---", file=sys.stderr)
    print(f"Contig: {contig_fname}", file=sys.stderr)
    print(f"Sylph taxon string: {taxon_string}", file=sys.stderr)
    if gtdb_record:
        print(f"GTDB organism name: {gtdb_record.get('ncbi_organism_name')}", file=sys.stderr)
        print(f"GTDB taxonomy: {gtdb_record.get('gtdb_taxonomy')}", file=sys.stderr)

    if not suggestions:
        print("No literal name match found.", file=sys.stderr)

    for i, (tax_id, name) in enumerate(suggestions, start=1):
        print(f"  [{i}] tax_id={tax_id} name={name}", file=sys.stderr)

    print("Enter a number to accept, 's' to skip, or 'q' to quit: ", file=sys.stderr, end="")
    sys.stderr.flush()

    while True:
        choice = tty.readline().strip()
        if choice.lower() == "q":
            print("Quitting at user request.", file=sys.stderr)
            sys.exit(1)
        if choice.lower() == "s" or choice == "":
            return None
        if choice.isdigit() and 1 <= int(choice) <= len(suggestions):
            return suggestions[int(choice) - 1]
        print("Invalid choice, try again: ", file=sys.stderr, end="")
        sys.stderr.flush()


def run(args):
    gtdb_lookup = load_gtdb_metadata(args.gtdb_metadata)
    parsed_taxonomy, name_alias_lookup, merged_map, deleted_ids = load_ncbi_taxonomy(
        args.taxonomy
    )

    report = ReportWriter(args.report)
    tty = Tty()

    tally = {
        "gtdb_accession": 0,
        "gtdb_species_taxid": 0,
        "literal_name": 0,
        "confirmed": 0,
        "unclassified": 0,
    }
    cache = {}

    try:
        with open(args.sylph_taxonomy, "r") as infile:
            reader = csv.DictReader(
                infile,
                delimiter="\t",
                fieldnames=["contig_fname", "taxon_string"],
            )

            for i, row in enumerate(reader):
                if i > 0 and i % 10000 == 0:
                    print(f"Processed {i} rows, current tally:\n{tally}", file=sys.stderr)

                contig_fname = row["contig_fname"]
                taxon_string = row["taxon_string"]

                cache_key = contig_fname
                if cache_key in cache:
                    tax_id, human_readable = cache[cache_key]
                    print(
                        "\t".join((contig_fname, taxon_string, tax_id, human_readable)),
                        file=sys.stdout,
                    )
                    continue

                gtdb_record = lookup_gtdb_record(gtdb_lookup, contig_fname)

                resolved_tax_id = None
                human_readable = None
                method = None

                if gtdb_record:
                    resolved_tax_id, status = resolve_tax_id(
                        gtdb_record.get("ncbi_taxid"), parsed_taxonomy, merged_map, deleted_ids
                    )
                    if resolved_tax_id:
                        human_readable = parsed_taxonomy[resolved_tax_id].human_readable
                        method = "gtdb_accession"
                        tally["gtdb_accession"] += 1
                    else:
                        species_tax_id, species_status = resolve_tax_id(
                            gtdb_record.get("ncbi_species_taxid"),
                            parsed_taxonomy,
                            merged_map,
                            deleted_ids,
                        )
                        if species_tax_id:
                            resolved_tax_id = species_tax_id
                            human_readable = parsed_taxonomy[species_tax_id].human_readable
                            method = "gtdb_species_taxid"
                            tally["gtdb_species_taxid"] += 1

                if not resolved_tax_id and gtdb_record:
                    matches = literal_name_match(
                        gtdb_record.get("ncbi_organism_name"), name_alias_lookup
                    )
                    if len(matches) == 1:
                        resolved_tax_id = matches[0]
                        human_readable = parsed_taxonomy[resolved_tax_id].human_readable
                        method = "literal_name"
                        tally["literal_name"] += 1

                if not resolved_tax_id:
                    species_name = species_name_from_taxon_string(taxon_string)
                    matches = literal_name_match(species_name, name_alias_lookup)
                    if len(matches) == 1:
                        resolved_tax_id = matches[0]
                        human_readable = parsed_taxonomy[resolved_tax_id].human_readable
                        method = "literal_name_taxon_string"
                        tally["literal_name"] += 1

                if resolved_tax_id:
                    cache[cache_key] = (resolved_tax_id, human_readable)
                    print(
                        "\t".join(
                            (contig_fname, taxon_string, resolved_tax_id, human_readable)
                        ),
                        file=sys.stdout,
                    )
                    continue

                # Unresolved: build suggestions from any available literal alias matches
                suggestions = []
                candidate_names = []
                if gtdb_record and gtdb_record.get("ncbi_organism_name"):
                    candidate_names.append(gtdb_record["ncbi_organism_name"])
                species_name = species_name_from_taxon_string(taxon_string)
                if species_name:
                    candidate_names.append(species_name)

                seen_tax_ids = set()
                for name in candidate_names:
                    for tax_id in literal_name_match(name, name_alias_lookup):
                        if tax_id not in seen_tax_ids:
                            seen_tax_ids.add(tax_id)
                            suggestions.append(
                                (tax_id, parsed_taxonomy[tax_id].human_readable)
                            )

                accepted_choice = confirm_match(
                    tty, taxon_string, contig_fname, gtdb_record, suggestions
                )

                failure_reason = "no_gtdb_record" if not gtdb_record else "taxid_absent_from_dump"

                report.write(
                    {
                        "contig_fname": contig_fname,
                        "taxon_string": taxon_string,
                        "gtdb_accession": gtdb_record.get("accession") if gtdb_record else "",
                        "gtdb_ncbi_taxid": gtdb_record.get("ncbi_taxid") if gtdb_record else "",
                        "gtdb_ncbi_organism_name": (
                            gtdb_record.get("ncbi_organism_name") if gtdb_record else ""
                        ),
                        "gtdb_taxonomy": gtdb_record.get("gtdb_taxonomy") if gtdb_record else "",
                        "failure_reason": failure_reason,
                        "suggested_tax_id": accepted_choice[0] if accepted_choice else (
                            suggestions[0][0] if suggestions else ""
                        ),
                        "suggested_name": accepted_choice[1] if accepted_choice else (
                            suggestions[0][1] if suggestions else ""
                        ),
                        "accepted": bool(accepted_choice),
                        "resolution_method": "confirmed_literal_name" if accepted_choice else "unclassified",
                    }
                )

                if accepted_choice:
                    tax_id, human_readable = accepted_choice
                    cache[cache_key] = (tax_id, human_readable)
                    tally["confirmed"] += 1
                    print(
                        "\t".join((contig_fname, taxon_string, tax_id, human_readable)),
                        file=sys.stdout,
                    )
                else:
                    cache[cache_key] = ("0", "Unclassified")
                    tally["unclassified"] += 1
                    print(
                        "\t".join((contig_fname, taxon_string, "0", "Unclassified")),
                        file=sys.stdout,
                    )
    finally:
        report.close()
        tty.close()

    print(f"Final tally:\n{tally}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Rationalise sylph taxon strings to NCBI taxIDs using GTDB metadata accession joins, "
        "falling back to literal (non-fuzzy) name matching with interactive confirmation."
    )
    parser.add_argument(
        "taxonomy",
        type=Path,
        help="Directory containing the NCBI taxonomy files (names.dmp, nodes.dmp, optionally merged.dmp/delnodes.dmp)",
    )
    parser.add_argument(
        "sylph_taxonomy",
        type=Path,
        help="Path to the sylph taxonomy tsv file (contig_fname, taxon_string columns)",
    )
    parser.add_argument(
        "--gtdb_metadata",
        type=Path,
        action="append",
        required=True,
        help="Path to a GTDB metadata TSV (e.g. bac120_metadata_r232.tsv.gz). May be given multiple times.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
        help="Path to write a TSV report of all taxa that could not be resolved via a direct GTDB/taxID match.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
