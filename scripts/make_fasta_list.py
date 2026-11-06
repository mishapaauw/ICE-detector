#!/usr/bin/env python
"""Build a PPanGGOLiN fasta list from a folder of genome assemblies.

Output is tab-separated: genome_name, absolute_path[, circular contig names...].
Genome name = file name without extension. Fails if a contig name occurs in
more than one genome, because PPanGGOLiN needs unique contig names.
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

FASTA_EXT = {".fasta", ".fa", ".fna", ".fas"}


def contig_names(path):
    with open(path) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]


def build_list(input_dir, circular=True):
    genomes = sorted(p for p in Path(input_dir).iterdir() if p.suffix in FASTA_EXT)
    if not genomes:
        sys.exit(f"No FASTA files ({', '.join(sorted(FASTA_EXT))}) found in {input_dir}")

    rows, seen = [], defaultdict(list)
    for path in genomes:
        contigs = contig_names(path)
        for c in contigs:
            seen[c].append(path.stem)
        row = [path.stem, str(path.resolve())]
        if circular:
            row += contigs
        rows.append(row)

    dups = {c: g for c, g in seen.items() if len(g) > 1}
    if dups:
        lines = "\n".join(f"  {c}: {', '.join(g)}" for c, g in sorted(dups.items()))
        sys.exit(f"Contig names must be unique across genomes. Duplicates:\n{lines}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_dir", help="folder with one FASTA per genome")
    ap.add_argument("-o", "--output", required=True, help="output fasta list (TSV)")
    ap.add_argument("--linear", action="store_true",
                    help="do not mark contigs as circular (default: all circular)")
    args = ap.parse_args()

    rows = build_list(args.input_dir, circular=not args.linear)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")
    print(f"Wrote {len(rows)} genomes to {args.output}")


if __name__ == "__main__":
    main()
