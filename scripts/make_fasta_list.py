#!/usr/bin/env python
"""Build a PPanGGOLiN fasta list from a folder of genome assemblies.

Output is tab-separated: genome_name, absolute_path[, circular contig names...].
Genome name = file name without extension. Fails if a contig name occurs in
more than one genome, because PPanGGOLiN needs unique contig names.

With --chromosome-dir, only the largest contig of each genome (the chromosome,
in a complete assembly) is used: it is written to <chromosome-dir>/<genome>.fasta
and the list points to these copies. The input files are only read. The copies
are temporary (run_all.py deletes them after the last step that needs them);
--report records which contig was kept and which were dropped.
"""
import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path

FASTA_EXT = {".fasta", ".fa", ".fna", ".fas"}


def contig_names(path):
    with open(path) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]


def read_fasta(path):
    """[(name, header line, sequence lines)] in file order."""
    records = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                records.append((line[1:].split()[0], line, []))
            elif records:
                records[-1][2].append(line)
    return records


def seq_length(lines):
    return sum(len(l.strip()) for l in lines)


def extract_chromosomes(input_dir, outdir):
    """Write the largest contig of each genome to outdir/<genome>.fasta.

    Ties: the first of the largest contigs in the file. Returns one report row
    per genome. Warns when the kept contig is < 50% of the assembly (draft genome?).
    """
    genomes = sorted(p for p in Path(input_dir).iterdir() if p.suffix in FASTA_EXT)
    if not genomes:
        sys.exit(f"No FASTA files ({', '.join(sorted(FASTA_EXT))}) found in {input_dir}")
    outdir = Path(outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    report = []
    for path in genomes:
        out = outdir / f"{path.stem}.fasta"
        if out.exists():
            sys.exit(f"Two input files give genome name {path.stem}")
        records = [(name, header, lines, seq_length(lines)) for name, header, lines in read_fasta(path)]
        if not records:
            sys.exit(f"No sequences in {path}")
        name, header, lines, length = max(records, key=lambda r: r[3])  # first max wins
        out.write_text(header + "".join(lines))
        total = sum(r[3] for r in records)
        dropped = [f"{n}:{l}" for n, _, _, l in records if n != name]
        if length < 0.5 * total:
            print(f"WARNING: {path.stem}: largest contig {name} is only {length / total:.0%} "
                  f"of the assembly; is this a complete genome?", file=sys.stderr)
        report.append({"genome": path.stem, "input_file": str(path.resolve()),
                       "kept_contig": name, "kept_length": length,
                       "assembly_length": total, "n_contigs": len(records),
                       "dropped_contigs": ",".join(dropped) or "-"})
    return report


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
    ap.add_argument("--chromosome-dir",
                    help="keep only the largest contig per genome, written to this folder")
    ap.add_argument("--report", help="TSV of kept/dropped contigs (with --chromosome-dir)")
    args = ap.parse_args()

    genome_dir = args.input_dir
    if args.chromosome_dir:
        report = extract_chromosomes(args.input_dir, args.chromosome_dir)
        genome_dir = args.chromosome_dir
        if args.report:
            with open(args.report, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(report[0]), delimiter="\t")
                w.writeheader()
                w.writerows(report)
        n_dropped = sum(r["n_contigs"] - 1 for r in report)
        print(f"Kept the largest contig of each genome; dropped {n_dropped} other contigs")
    rows = build_list(genome_dir, circular=not args.linear)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")
    print(f"Wrote {len(rows)} genomes to {args.output}")


if __name__ == "__main__":
    main()
