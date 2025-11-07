#!/usr/bin/env python3
import sys
import argparse
import re


def get_gene_name(attr_field: str) -> str:
    """Extract gene_name or gene_id from GTF attributes using regex"""
    # Try gene_name first
    m = re.search(r'gene_name "([^"]+)"', attr_field)
    if m:
        return m.group(1)
    # Fallback to gene_id
    m = re.search(r'gene_id "([^"]+)"', attr_field)
    if m:
        return m.group(1)
    return "NA"


def gtf_to_bed(gtf_path: str, bed_out: str):
    with open(gtf_path, "r") as fin, open(bed_out, "w") as fout:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "gene":
                continue
            try:
                start0 = max(0, int(start) - 1)  # GTF is 1-based inclusive; BED is 0-based half-open
                end0 = int(end)
            except ValueError:
                continue
            gene = get_gene_name(attrs)
            fout.write(f"{chrom}\t{start0}\t{end0}\t{gene}\t{strand}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert GTF (gene features) to BED: chrom,start,end,gene,strand")
    parser.add_argument("gtf", help="Input GTF file (gencode.v49.basic.annotation.gtf)")
    parser.add_argument("bed_out", help="Output BED file path")
    args = parser.parse_args()

    gtf_to_bed(args.gtf, args.bed_out)


if __name__ == "__main__":
    main()


