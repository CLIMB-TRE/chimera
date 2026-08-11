#!/usr/bin/env python3

import argparse
import os
import re
import sys
import tempfile
import pysam

XA_RE = re.compile(r"([^,]+),([+-])(\d+),([^,]+),(\d+)")
CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")
QUERY_CONSUMING_OPS = set("MIS=X")


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def cigar_query_length(cigar: str) -> int:
    """Sum of CIGAR operations that consume query bases (M/I/S/=/X), i.e. the
    read length implied by the CIGAR, excluding hard-clipped (H) bases."""
    return sum(
        int(length) for length, op in CIGAR_OP_RE.findall(cigar) if op in QUERY_CONSUMING_OPS
    )


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

    if primary.query_sequence is None or cigar_query_length(cigar) != len(
        primary.query_sequence
    ):
        print(
            f"Skipping XA hit for {primary.query_name} on {ref_name}: CIGAR {cigar} "
            f"implies a read length that does not match the primary alignment's sequence.",
            file=sys.stderr,
        )
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
    if primary.query_sequence is None or not record.cigarstring or cigar_query_length(
        record.cigarstring
    ) != len(primary.query_sequence):
        print(
            f"Skipping SEQ fill for secondary alignment of {record.query_name}: CIGAR "
            f"{record.cigarstring} does not match the primary alignment's sequence length.",
            file=sys.stderr,
        )
        return record

    if record.is_reverse != primary.is_reverse:
        record.query_sequence = reverse_complement(primary.query_sequence)
        if primary.query_qualities is not None:
            record.query_qualities = primary.query_qualities[::-1]
    else:
        record.query_sequence = primary.query_sequence
        record.query_qualities = primary.query_qualities
    return record


def process_group(records: list, header: pysam.AlignmentHeader, out_bam: pysam.AlignmentFile):
    """Process every record for a single query name (all mates, primary and
    secondary/supplementary alike) together, then write them all out."""
    primaries = {}
    others = []
    for read in records:
        if not read.is_secondary and not read.is_supplementary:
            primaries[(read.is_read1, read.is_read2)] = read
        else:
            others.append(read)

    for read in primaries.values():
        out_bam.write(read)
        if read.has_tag("XA"):
            xa = read.get_tag("XA")
            if xa:
                for ref_name, strand, pos, cigar, nm in parse_xa_tag(xa):
                    secondary = make_secondary_from_xa(read, header, ref_name, strand, pos, cigar, nm)
                    if secondary is not None:
                        out_bam.write(secondary)

    for record in others:
        primary = primaries.get((record.is_read1, record.is_read2))
        if primary is not None and (
            record.query_sequence is None or record.query_sequence == ""
        ):
            record = fill_secondary_seq(record, primary)
        out_bam.write(record)


def run(args):
    # Name-collate first so that every record for a read (primary, secondary,
    # supplementary; both mates) ends up adjacent. samtools collate buckets
    # records to temporary files rather than sorting in memory, so this scales
    # to large inputs; it also means we only ever need to hold one query
    # name's worth of records in memory below, instead of the whole file.
    with tempfile.TemporaryDirectory(dir=".") as tmp_dir:
        collated_bam = os.path.join(tmp_dir, "collated.bam")
        pysam.collate("--no-PG", "-o", collated_bam, str(args.input_bam))

        with pysam.AlignmentFile(collated_bam, "rb") as in_bam:
            header = in_bam.header
            with pysam.AlignmentFile(args.output_bam, "wb", header=header) as out_bam:
                current_qname = None
                group = []
                for read in in_bam:
                    if read.query_name != current_qname:
                        if group:
                            process_group(group, header, out_bam)
                        group = []
                        current_qname = read.query_name
                    group.append(read)
                if group:
                    process_group(group, header, out_bam)


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
