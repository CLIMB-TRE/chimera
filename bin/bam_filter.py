#!/usr/bin/env python3

import pysam

import itertools
from collections import defaultdict


def read_pair_generator(reads, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads (including secondary alignments) are queued by (query_name, is_secondary)
    and paired off in encounter order, so that every secondary alignment is matched
    with a mate where one exists. Any record left unmatched once the input is
    exhausted (e.g. a secondary present for only one mate) is yielded on its own
    rather than being silently dropped.
    """
    read1_queues = defaultdict(list)
    read2_queues = defaultdict(list)
    for read in reads:
        if read.is_supplementary:
            continue
        key = (read.query_name, read.is_secondary)
        if read.is_read1:
            read1_queues[key].append(read)
        else:
            read2_queues[key].append(read)

        while read1_queues[key] and read2_queues[key]:
            yield read1_queues[key].pop(0), read2_queues[key].pop(0)

    for key, queue in read1_queues.items():
        for read in queue:
            yield read, None
    for key, queue in read2_queues.items():
        for read in queue:
            yield None, read


def read_checks(
    read: pysam.AlignedSegment | tuple[pysam.AlignedSegment, pysam.AlignedSegment],
    min_alignment_proportion: float,
    verbose: bool = False,
) -> bool:
    if isinstance(read, tuple):
        read1, read2 = read

        if read1 is not None and read2 is not None and read1.query_name != read2.query_name:
            return False

        if read1 is None:
            return read_checks(read2, min_alignment_proportion, verbose)
        if read2 is None:
            return read_checks(read1, min_alignment_proportion, verbose)

        if read1.is_unmapped or read2.is_unmapped:
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
        default=0.0,
        help="Minimum alignment proportion (how much of the source read is aligned) to retain a read (0-1)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    with (
        pysam.AlignmentFile(args.input_bam, "rb") as in_bam,
        pysam.AlignmentFile(args.output_bam, "wb", template=in_bam) as out_bam,
    ):

        try:
            first_read = next(in_bam)
        except StopIteration:
            return

        chained_iterator = itertools.chain([first_read], in_bam)

        if first_read.is_paired:
            read_pairs = read_pair_generator(chained_iterator)

        for read in read_pairs if first_read.is_paired else chained_iterator:
            if read_checks(read, args.min_alignment_proportion, args.verbose):
                if first_read.is_paired:
                    for mate in read:
                        if mate is not None:
                            out_bam.write(mate)
                else:
                    out_bam.write(read)


if __name__ == "__main__":
    main()
