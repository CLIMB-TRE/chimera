#!/usr/bin/env python3

import numpy as np
import math
import csv
import sys
import pysam


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
                "start_end_positions": {},
                "num_reads": 0,
                "forward_reads": 0,
                "unique_mappers": 0,
                "primary_alignments": 0,
            },
        )
        stats_dict[ref_name]["num_reads"] += 1

        try:
            nm_tag = int(read.get_tag("NM"))
        except Exception:
            print("NM tag not found for read, exiting:")
            print(read)
            sys.exit(1)

        read_ref_map.setdefault(read.query_name, set())
        read_ref_map[read.query_name].add(ref_name)

        aln_length = read.query_alignment_length
        identity = ((aln_length - nm_tag) / aln_length) * 100

        start_end_tuple = (read.reference_start, read.reference_end)
        stats_dict[ref_name]["start_end_positions"].setdefault(start_end_tuple, 0)
        stats_dict[ref_name]["start_end_positions"][start_end_tuple] += 1

        if not read.is_secondary:
            stats_dict[ref_name]["primary_alignments"] += 1

        stats_dict[ref_name]["identities"].append(identity)
        stats_dict[ref_name]["alignment_lengths"].append(aln_length)
        stats_dict[ref_name]["read_lengths"].append(read.query_length)
        stats_dict[ref_name]["alignment_proportions"].append(
            aln_length / read.query_length
        )

        if not read.is_reverse:
            stats_dict[ref_name]["forward_reads"] += 1

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
            "mean_identity": mean_identity if mean_identity > 0 else 0,
            "duplication_rate": duplication_rate if duplication_rate > 0 else 0,
            "mean_aln_length": mean_aln_length if mean_aln_length > 0 else 0,
            "forward_proportion": forward_proportion if forward_proportion > 0 else 0,
            "uniquely_mapped_reads": stats["unique_mappers"],
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


def coverage_tsv_parser(depth_tsv: str) -> dict:
    """
    Parse a depth TSV file and return a dictionary with coverage information.

    Parameters
    ----------
    depth_tsv : str
        Path to the depth TSV file.

    Returns
    -------
    list
        A list of dictionaries, each containing coverage information for a specific region.
    """
    coverage_dict = {}

    with open(depth_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            coverage_dict[row["#rname"]] = row

    return coverage_dict


def alignment_stats(depth_array: np.ndarray, coverage_stats: dict) -> dict:
    """Generate some basic stats from a depth array, including coverage evenness (E), mean depth, breadth at 1x and 10x, mapped reads, mapped bases.

    Args:
        depth_array (np.ndarray): Array of per-base coverage values.

    Returns:
        dict: A dictionary containing the computed alignment statistics.
    """
    stats = {
        "evenness_value": coverage_evenness(depth_array),
        "mean_depth": int(depth_array.mean()),
        "coverage_1x": int((depth_array > 0).sum() / len(depth_array) * 100),
        "coverage_10x": int((depth_array > 9).sum() / len(depth_array) * 100),
        "mapped_reads": int(coverage_stats["numreads"]),
        "mapped_bases": int(
            int(coverage_stats["endpos"]) * float(coverage_stats["meandepth"])
        ),
        "mean_read_length": (
            int(np.mean(coverage_stats["read_lengths"]))
            if coverage_stats["read_lengths"]
            else 0
        ),
        "mean_alignment_proportion": (
            float(np.mean(coverage_stats["alignment_proportions"]))
            if coverage_stats["alignment_proportions"]
            else 0.0
        ),
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


def run(args):

    depth_arrays = depth_tsv_to_np_arrays(args.depth_tsv)
    coverage_info = coverage_tsv_parser(args.coverage_tsv)
    reference_metadata = reference_metadata_parser(args.database_metadata)
    bam_stats = generate_bam_stats(args.bam)

    # Print report
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
            "primary_alignments",
            "mapped_bases",
            "mean_read_identity",
            "read_duplication_rate",
            "forward_proportion",
            "mean_read_length",
            "mean_alignment_length",
            "mean_alignment_proportion",
        ],
    )
    writer.writeheader()

    for ref in depth_arrays:
        if ref in coverage_info:
            stats = alignment_stats(depth_arrays[ref], coverage_info[ref])
            stats["unique_accession"] = ref
            stats["taxon_id"] = reference_metadata[ref]["taxon_id"]
            stats["human_readable"] = reference_metadata[ref]["human_readable"]
            stats["accession_description"] = reference_metadata[ref][
                "accession_description"
            ]
            stats["sequence_length"] = reference_metadata[ref]["sequence_length"]
        else:
            print(f"ERROR: Reference {ref} found in depth TSV but not in coverage TSV.")
            sys.exit(1)

        if ref in bam_stats:
            stats["mean_read_identity"] = bam_stats[ref]["mean_identity"]
            stats["read_duplication_rate"] = bam_stats[ref]["duplication_rate"]
            stats["mean_alignment_length"] = bam_stats[ref]["mean_aln_length"]
            stats["forward_proportion"] = bam_stats[ref]["forward_proportion"]
            stats["uniquely_mapped_reads"] = bam_stats[ref]["uniquely_mapped_reads"]
            stats["primary_alignments"] = bam_stats[ref]["primary_alignments"]
            stats["mean_read_length"] = bam_stats[ref]["mean_read_length"]
            stats["mean_alignment_proportion"] = bam_stats[ref][
                "mean_alignment_proportion"
            ]
        else:
            print(f"WARNING: Reference {ref} found in depth TSV but not in BAM stats.")
            sys.exit(1)

        writer.writerow(stats)


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
        "--coverage_tsv",
        type=str,
        required=True,
        help="Path to the coverage TSV file, generated by samtools coverage.",
    )
    parser.add_argument(
        "--database_metadata",
        type=str,
        required=True,
        help="Path to the database metadata TSV file, containing reference taxonomy etc.",
    )
    parser.add_argument("--bam", type=str, required=True, help="Path to the BAM file.")
    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
