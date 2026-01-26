#!/usr/bin/env python3

import pysam

import itertools
from collections import defaultdict


def read_pair_generator(reads, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in reads:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_checks(
    read: pysam.AlignedSegment | tuple[pysam.AlignedSegment, pysam.AlignedSegment],
    min_alignment_proportion: float,
    verbose: bool = False,
) -> bool:
    if isinstance(read, tuple):
        read1, read2 = read

        if read1 is None or read2 is None:
            return False

        if read1.is_unmapped or read2.is_unmapped:
            return False

        if read1.query_name != read2.query_name:
            return False

        aligned_length1 = read1.query_alignment_length
        total_length1 = read1.infer_read_length()

        aligned_length2 = read2.query_alignment_length
        total_length2 = read2.infer_read_length()

        prop1 = 0 if total_length1 == 0 else aligned_length1 / total_length1
        prop2 = 0 if total_length2 == 0 else aligned_length2 / total_length2

        if prop1 < min_alignment_proportion and prop2 < min_alignment_proportion:
            return False

    else:
        if read.is_unmapped:
            return False

        aligned_length = read.query_alignment_length
        total_length = read.infer_read_length()

        if (
            total_length == 0
            or (aligned_length / total_length) < min_alignment_proportion
        ):
            return False
        else:
            if verbose:
                print(
                    f"Keeping read {read.query_name} aligned to {read.reference_name} with alignment proportion {aligned_length / total_length}"
                )
                print(f"Aligned length: {aligned_length}, Total length: {total_length}")

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Filter BAM file by mapping quality.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument(
        "--min_alignment_proportion",
        type=float,
        help="Minimum alignment proportion (how much of the source read is aligned) to retain a read (0-1)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    with (
        pysam.AlignmentFile(args.input_bam, "rb") as in_bam,
        pysam.AlignmentFile(args.output_bam, "wb", template=in_bam) as out_bam,
    ):

        first_read = next(in_bam)
        chained_iterator = itertools.chain([first_read], in_bam)

        if first_read.is_paired:
            read_pairs = read_pair_generator(chained_iterator)

        for read in read_pairs if first_read.is_paired else chained_iterator:
            if read_checks(read, args.min_alignment_proportion, args.verbose):
                if first_read.is_paired:
                    out_bam.write(read[0])
                    out_bam.write(read[1])
                else:
                    out_bam.write(read)


if __name__ == "__main__":
    main()
