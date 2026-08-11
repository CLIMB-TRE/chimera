#!/usr/bin/env python3

import argparse
import re
import pysam

XA_RE = re.compile(r"([^,]+),([+-])(\d+),([^,]+),(\d+)")


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def parse_xa_tag(xa: str) -> list[tuple[str, str, int, str, int]]:
    """Parse a bwa-mem2 XA tag into (ref_name, strand, pos, cigar, nm) tuples."""
    hits = []
    for entry in xa.rstrip(";").split(";"):
        if not entry:
            continue
        match = XA_RE.fullmatch(entry)
        if not match:
            continue
        ref_name, strand, pos, cigar, nm = match.groups()
        hits.append((ref_name, strand, int(pos), cigar, int(nm)))
    return hits


def make_secondary_from_xa(
    primary: pysam.AlignedSegment,
    header: pysam.AlignmentHeader,
    ref_name: str,
    strand: str,
    pos: int,
    cigar: str,
    nm: int,
) -> pysam.AlignedSegment | None:
    reference_id = header.get_tid(ref_name)
    if reference_id == -1:
        return None

    secondary = pysam.AlignedSegment(header)
    secondary.query_name = primary.query_name
    secondary.reference_id = reference_id
    secondary.reference_start = pos - 1
    secondary.mapping_quality = 0
    secondary.cigarstring = cigar

    flag = primary.flag & ~0x100
    flag |= 0x100
    is_reverse = strand == "-"
    if is_reverse != primary.is_reverse:
        flag ^= 0x10

    secondary.flag = flag

    if is_reverse != primary.is_reverse:
        secondary.query_sequence = reverse_complement(primary.query_sequence)
        secondary.query_qualities = primary.query_qualities[::-1]
    else:
        secondary.query_sequence = primary.query_sequence
        secondary.query_qualities = primary.query_qualities

    secondary.set_tag("NM", nm)
    secondary.next_reference_id = primary.next_reference_id
    secondary.next_reference_start = primary.next_reference_start
    secondary.template_length = 0

    return secondary


def fill_secondary_seq(record: pysam.AlignedSegment, primary: pysam.AlignedSegment):
    """Fill SEQ/QUAL on an existing secondary record from its own primary record."""
    if record.is_reverse != primary.is_reverse:
        record.query_sequence = reverse_complement(primary.query_sequence)
        if primary.query_qualities is not None:
            record.query_qualities = primary.query_qualities[::-1]
    else:
        record.query_sequence = primary.query_sequence
        record.query_qualities = primary.query_qualities
    return record


def run(args):
    with pysam.AlignmentFile(args.input_bam, "rb") as in_bam:
        header = in_bam.header
        with pysam.AlignmentFile(args.output_bam, "wb", header=header) as out_bam:
            primaries_by_read = {}
            secondaries_by_read = {}

            for read in in_bam:
                key = (read.query_name, read.is_read1, read.is_read2)
                if not read.is_secondary and not read.is_supplementary:
                    primaries_by_read[key] = read
                    out_bam.write(read)
                    xa = None
                    if read.has_tag("XA"):
                        xa = read.get_tag("XA")
                    if xa:
                        for ref_name, strand, pos, cigar, nm in parse_xa_tag(xa):
                            secondary = make_secondary_from_xa(
                                read, header, ref_name, strand, pos, cigar, nm
                            )
                            if secondary is not None:
                                out_bam.write(secondary)
                else:
                    secondaries_by_read.setdefault(key, []).append(read)

            for key, records in secondaries_by_read.items():
                primary = primaries_by_read.get(key)
                for record in records:
                    if primary is not None and (
                        record.query_sequence is None or record.query_sequence == ""
                    ):
                        record = fill_secondary_seq(record, primary)
                    out_bam.write(record)


def main():
    parser = argparse.ArgumentParser(
        description="Expand bwa-mem2 XA tags into secondary alignment records and fill SEQ/QUAL on any secondary records missing sequence, using the primary alignment of the same read as the source."
    )
    parser.add_argument("input_bam", help="Input BAM file (from bwa-mem2)")
    parser.add_argument("output_bam", help="Output BAM file")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
