#!/usr/bin/env python3

import numpy as np
import math
import csv
import sys


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
        "1x_coverage": int((depth_array > 0).sum() / len(depth_array) * 100),
        "10x_coverage": int((depth_array > 9).sum() / len(depth_array) * 100),
        "mapped_reads": int(coverage_stats["numreads"]),
        "mapped_bases": int(
            int(coverage_stats["endpos"]) * float(coverage_stats["meandepth"])
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

    out_rows = 

    for ref in depth_arrays:
        if ref in coverage_info:
            stats = alignment_stats(depth_arrays[ref], coverage_info[ref])
            stats["reference"] = ref
            out_rows.append(stats)
        else:
            print(f"ERROR: Reference {ref} found in depth TSV but not in coverage TSV.")
            sys.exit(1)

    # Print report
    writer = csv.DictWriter(
        sys.stdout,
        delimiter="\t",
        fieldnames=[
            "reference",
            "tax_id",
            "human_readable",            
            "evenness_value",
            "mean_depth",
            "1x_coverage",
            "10x_coverage",
            "mapped_reads",
            "mapped_bases",
        ],
    )
    writer.writeheader()

    for row in out_rows:
        row["tax_id"] = reference_metadata[row["reference"]]["tax_id"]
        row["human_readable"] = reference_metadata[row["reference"]]["human_readable"]
        writer.writerow(row)


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
    args = parser.parse_args()

    run(args)
