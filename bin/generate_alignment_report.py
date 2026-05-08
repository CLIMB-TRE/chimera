#!/usr/bin/env python3

import numpy as np
import math
import csv
import sys
import pysam
import json
import jsonschema


def generate_alignment_complexity(read: pysam.AlignedSegment) -> float:
    """
    Calculate the complexity of an aligned section of a read based on the proportion of consecutive identical bases.
    e.g. How many bases are the same as the previous base, divided by total bases - 1. AAAGAGA would have a complexity score of 0.285... (2/7).

    Parameters
    ----------
    read : pysam.AlignedSegment
        A pysam AlignedSegment object representing a read.

    Returns
    -------
    float
        Complexity score between 0 and 1, where 1 indicates high complexity and 0 indicates low complexity.
    """

    seq = read.query_alignment_sequence
    if not seq or len(seq) < 2:
        return 0.0  # No sequence or too short to determine complexity

    same_base_count = sum(1 for i in range(1, len(seq)) if seq[i] == seq[i - 1])
    complexity = 1 - (same_base_count / (len(seq) - 1))

    return complexity


def generate_bam_stats(bam_file: str) -> dict:
    """
    Get basic stats for each CHROM from a BAM file, specifically:
       1. the mean read BLAST-like identity with reference sequences.
       2. Duplication rate (percentage of reads which start and end at the same position as another read).
       3. Mean alignment length.
       4. Forward strand proportion (percentage of reads mapped to the forward strand) -> strand bias.

    Parameters
    ----------
    bam_file : str
        Path to the BAM file.

    Returns
    -------
    dict
        A dictionary where keys are reference names and values are stats dictionaries.
    """
    stats_dict = {}
    read_ref_map = {}

    bam = pysam.AlignmentFile(bam_file, "rb")

    for read in bam:
        if read.is_unmapped:
            continue

        ref_name = bam.get_reference_name(read.reference_id)

        stats_dict.setdefault(
            ref_name,
            {
                "identities": [],
                "alignment_lengths": [],
                "read_lengths": [],
                "alignment_proportions": [],
                "alignment_complexities": [],
                "start_end_positions": {},
                "num_reads": 0,
                "forward_reads": 0,
                "unique_mappers": 0,
            },
        )
        stats_dict[ref_name]["num_reads"] += 1

        try:
            nm_tag = int(read.get_tag("NM"))
        except Exception:
            print("NM tag not found for read, exiting:")
            print(read)
            sys.exit(1)

        try:
            read_ref_map.setdefault(read.query_name, set())
            read_ref_map[read.query_name].add(ref_name)

            aln_length = read.query_alignment_length
            identity = ((aln_length - nm_tag) / aln_length) * 100

            start_end_tuple = (read.reference_start, read.reference_end)
            stats_dict[ref_name]["start_end_positions"].setdefault(start_end_tuple, 0)
            stats_dict[ref_name]["start_end_positions"][start_end_tuple] += 1

            stats_dict[ref_name]["identities"].append(identity)
            stats_dict[ref_name]["alignment_lengths"].append(aln_length)
            stats_dict[ref_name]["read_lengths"].append(read.infer_read_length())
            stats_dict[ref_name]["alignment_proportions"].append(
                aln_length / read.infer_read_length()
            )
            stats_dict[ref_name]["alignment_complexities"].append(
                generate_alignment_complexity(read)
            )

            if not read.is_reverse:
                stats_dict[ref_name]["forward_reads"] += 1

        except Exception as e:
            print(f"Error processing read:\n{read}\nError: {e}", file=sys.stderr)
            sys.exit(1)

    bam.close()

    out_stats = {}

    for read_name in read_ref_map:
        if len(read_ref_map[read_name]) == 1:
            unique_ref = next(iter(read_ref_map[read_name]))
            stats_dict[unique_ref]["unique_mappers"] += 1

    for ref in stats_dict:
        stats = stats_dict[ref]
        duplicates = sum(x for x in stats["start_end_positions"].values() if x > 1)
        if duplicates > 0:
            duplication_rate = round(
                (duplicates / stats["num_reads"] * 100),
                2,
            )
        else:
            duplication_rate = 0

        mean_identity = round(np.mean(stats["identities"]), 2)
        mean_aln_length = round(np.mean(stats["alignment_lengths"]), 2)
        forward_proportion = round(
            stats["forward_reads"] / stats["num_reads"],
            2,
        )

        out_stats[ref] = {
            "num_reads": stats["num_reads"],
            "mean_identity": mean_identity if mean_identity > 0 else 0,
            "duplication_rate": duplication_rate if duplication_rate > 0 else 0,
            "mean_aln_length": mean_aln_length if mean_aln_length > 0 else 0,
            "forward_proportion": forward_proportion if forward_proportion > 0 else 0,
            "uniquely_mapped_reads": stats["unique_mappers"],
            "mean_read_length": round(np.mean(stats["read_lengths"]), 2),
            "mean_alignment_proportion": round(
                np.mean(stats["alignment_proportions"]), 2
            ),
            "mean_alignment_complexity": round(
                np.mean(stats["alignment_complexities"]), 2
            ),
        }

    return out_stats


def coverage_evenness(coverage: np.ndarray, num_points: int = 1000) -> int:
    """
    Compute coverage evenness score E (%) from per-base coverage.

    Parameters
    ----------
    coverage : numpy array
        Per-base coverage values (integers or floats).
    num_points : int, optional
        Number of grid points for numerical integration (default=1000).

    Returns
    -------
    E : int
        Evenness score as a percentage (0-100). Higher values indicate more even coverage.
        0% means all coverage is concentrated in a single position, while 100% means coverage is evenly distributed across all positions.
    """

    mean_cov = coverage.mean()

    if mean_cov == 0:
        return 0  # no coverage at all

    # Normalize coverage by mean
    norm_cov = coverage / mean_cov

    # Build thresholds between 0 and 1
    i_values = np.linspace(0, 1, num_points)

    # For each threshold i, compute fraction of bases >= i
    F_values = np.array([(norm_cov >= i).mean() for i in i_values])

    # Integrate F(i) over [0,1] using trapezoidal rule
    integral = np.trapezoid(F_values, i_values)

    # Convert to percentage
    E = math.floor(integral * 100)

    return E


def depth_tsv_to_np_arrays(depth_tsv: str) -> dict:
    """
    Read a depth TSV file and return positions and coverage

    Parameters
    ----------
    depth_tsv : str
        Path to the depth TSV file. The file should have three columns: reference name, position (1-based), and coverage. As produced by `samtools depth` with `-a` so that 0 depth positions are included.

    Returns
    -------
    dict
        A dictionary where keys are reference names and values 1d numpy arrays of coverage values.
    """

    depth_dict = {}
    depth_arrays = {}

    with open(depth_tsv, "r") as f:
        for line in f:
            ref, pos, cov = line.strip().split("\t")
            pos = int(pos) - 1  # Convert to 0-based index
            cov = int(cov)

            depth_dict.setdefault(ref, [])

            # Extend the list to the current position if necessary
            while len(depth_dict[ref]) <= pos:
                depth_dict[ref].append(0.0)

            depth_dict[ref][pos] = cov

    # Convert lists to numpy arrays
    for ref in depth_dict:
        depth_arrays[ref] = np.array(depth_dict[ref])

    return depth_arrays



def alignment_stats(depth_array: np.ndarray, num_reads: int) -> dict:
    """Generate some basic stats from a depth array, including coverage evenness (E), mean depth, breadth at 1x and 10x, mapped reads, mapped bases.

    Args:
        depth_array (np.ndarray): Array of per-base coverage values.
        num_reads (int): Number of mapped reads for this reference.

    Returns:
        dict: A dictionary containing the computed alignment statistics.
    """
    stats = {
        "evenness_value": coverage_evenness(depth_array),
        "mean_depth": int(depth_array.mean()),
        "coverage_1x": int((depth_array > 0).sum() / len(depth_array) * 100),
        "coverage_10x": int((depth_array > 9).sum() / len(depth_array) * 100),
        "mapped_reads": num_reads,
        "mapped_bases": int(depth_array.sum()),
    }

    return stats


def reference_metadata_parser(database_metadata: str) -> dict:
    """
    Parse a database metadata TSV file and return a dictionary with reference information.

    Parameters
    ----------
    database_metadata : str
        Path to the database metadata TSV file.

    Returns
    -------
    dict
        A dictionary where keys are reference names and values are metadata dictionaries.
    """
    metadata_dict = {}

    with open(database_metadata, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metadata_dict[row["unique_accession"]] = row

    return metadata_dict


def validate_scoring_matrix(scoring_matrix_path: str, json_schema_path: str) -> dict:
    """Validate the scoring matrix against a JSON schema and check range continuity.

    Args:
        scoring_matrix_path (str): path to the scoring matrix file
        json_schema_path (str): path to the JSON schema file

    Raises:
        jsonschema.ValidationError: If the scoring matrix does not conform to the JSON schema
        jsonschema.SchemaError: If the JSON schema itself is invalid
        ValueError: If any range has min >= max
        ValueError: If ranges are not contiguous
        ValueError: If an open-ended range is not the last one
    Returns:
        dict: The validated scoring matrix
    """

    with open(json_schema_path, "r") as schema_file:
        schema = json.load(schema_file)

    with open(scoring_matrix_path, "r") as matrix_file:
        scoring_matrix = json.load(matrix_file)

    jsonschema.validate(instance=scoring_matrix, schema=schema)

    for metric, details in scoring_matrix["metrics"].items():

        for i, r in enumerate(details["bands"]):
            if r["max"] is not None and r["min"] >= r["max"]:
                raise ValueError(f"{metric}: min must be < max")

            if i > 0:
                prev = details["bands"][i - 1]
                if prev["max"] != r["min"]:
                    raise ValueError(f"{metric}: ranges must be contiguous")

            if r["max"] is None and i != len(details["bands"]) - 1:
                raise ValueError(f"{metric}: open-ended range must be last")

    for i, r in enumerate(scoring_matrix["total_score_thresholds"]):
        if r["max"] is not None and r["min"] >= r["max"]:
            raise ValueError(f"Score category {r['label']}: min must be < max")

        if i > 0:
            prev = scoring_matrix["total_score_thresholds"][i - 1]
            if prev["max"] != r["min"]:
                raise ValueError(
                    f"Score category {r['label']}: ranges must be contiguous"
                )

        if r["max"] is None and i != len(scoring_matrix["total_score_thresholds"]) - 1:
            raise ValueError(
                f"Score category {r['label']}: open-ended range must be last"
            )

    return scoring_matrix


def score_record(record: dict, scoring_matrix: dict) -> int:
    total_score = 0

    for metric, details in scoring_matrix["metrics"].items():
        value = record.get(metric)
        if value is None:
            continue

        for band in details["bands"]:
            if band["max"] is None:
                if float(value) >= band["min"]:
                    total_score += band["score"]
                    break
            else:
                if band["min"] <= float(value) < band["max"]:
                    total_score += band["score"]
                    break

    return total_score


def total_score_category(record: dict, scoring_matrix: dict) -> str:
    score = score_record(record, scoring_matrix)

    thresholds = scoring_matrix["total_score_thresholds"]
    for category in thresholds:
        if category["max"] is None:
            if float(score) >= category["min"]:
                return category["label"]
        else:
            if category["min"] <= float(score) < category["max"]:
                return category["label"]


def run(args):

    depth_arrays = depth_tsv_to_np_arrays(args.depth_tsv)
    reference_metadata = reference_metadata_parser(args.database_metadata)
    bam_stats = generate_bam_stats(args.bam)

    ref_stat_rows = []

    for ref in depth_arrays:
        if ref not in bam_stats:
            print(f"ERROR: Reference {ref} found in depth TSV but not in BAM stats.")
            sys.exit(1)

        stats = alignment_stats(depth_arrays[ref], bam_stats[ref]["num_reads"])
        stats["unique_accession"] = ref
        stats["taxon_id"] = reference_metadata[ref]["taxon_id"]
        stats["human_readable"] = reference_metadata[ref]["human_readable"]
        stats["accession_description"] = reference_metadata[ref]["accession_description"]
        stats["sequence_length"] = reference_metadata[ref]["sequence_length"]
        stats["mean_read_identity"] = bam_stats[ref]["mean_identity"]
        stats["read_duplication_rate"] = bam_stats[ref]["duplication_rate"]
        stats["mean_alignment_length"] = bam_stats[ref]["mean_aln_length"]
        stats["forward_proportion"] = bam_stats[ref]["forward_proportion"]
        stats["uniquely_mapped_reads"] = bam_stats[ref]["uniquely_mapped_reads"]
        stats["mean_read_length"] = bam_stats[ref]["mean_read_length"]
        stats["mean_alignment_proportion"] = bam_stats[ref]["mean_alignment_proportion"]
        stats["mean_alignment_complexity"] = bam_stats[ref]["mean_alignment_complexity"]
        ref_stat_rows.append(stats)

    if not args.scoring_matrix:
        writer = csv.DictWriter(
            sys.stdout,
            delimiter="\t",
            fieldnames=[
                "taxon_id",
                "human_readable",
                "unique_accession",
                "accession_description",
                "sequence_length",
                "evenness_value",
                "mean_depth",
                "coverage_1x",
                "coverage_10x",
                "mapped_reads",
                "uniquely_mapped_reads",
                "mapped_bases",
                "mean_read_identity",
                "read_duplication_rate",
                "forward_proportion",
                "mean_read_length",
                "mean_alignment_length",
                "mean_alignment_proportion",
                "mean_alignment_complexity",
            ],
        )
        writer.writeheader()
        writer.writerows(ref_stat_rows)
        return

    if args.json_schema:
        scoring_matrix = validate_scoring_matrix(args.scoring_matrix, args.json_schema)
    else:
        print(
            "WARNING: No JSON schema provided so skipping scoring matrix validation. This may break scoring!",
            file=sys.stderr,
        )
        scoring_matrix = json.load(open(args.scoring_matrix, "r"))

    writer = csv.DictWriter(
        sys.stdout,
        delimiter="\t",
        fieldnames=[
            "taxon_id",
            "human_readable",
            "unique_accession",
            "accession_description",
            "sequence_length",
            "evenness_value",
            "mean_depth",
            "coverage_1x",
            "coverage_10x",
            "mapped_reads",
            "uniquely_mapped_reads",
            "mapped_bases",
            "mean_read_identity",
            "read_duplication_rate",
            "forward_proportion",
            "mean_read_length",
            "mean_alignment_length",
            "mean_alignment_proportion",
            "mean_alignment_complexity",
            "total_score",
            "confidence",
        ],
    )
    writer.writeheader()

    for row in ref_stat_rows:
        total_score = score_record(row, scoring_matrix)
        score_category = total_score_category(row, scoring_matrix)

        row["total_score"] = total_score
        row["confidence"] = score_category

    sorted_rows = sorted(
        ref_stat_rows,
        key=lambda x: x["total_score"],
        reverse=True,
    )

    writer.writerows(sorted_rows)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate alignment report from depth TSV file."
    )
    parser.add_argument(
        "--depth_tsv",
        type=str,
        required=True,
        help="Path to the depth TSV file, generated by samtools depth -a.",
    )
    parser.add_argument(
        "--database_metadata",
        type=str,
        required=True,
        help="Path to the database metadata TSV file, containing reference taxonomy etc.",
    )
    parser.add_argument(
        "--scoring_matrix",
        type=str,
        help="Path to the scoring matrix file.",
    )
    parser.add_argument(
        "--json_schema",
        type=str,
        help="Path to the JSON schema file for validating the scoring matrix.",
    )

    parser.add_argument("bam", type=str, help="Path to the BAM file.")
    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
