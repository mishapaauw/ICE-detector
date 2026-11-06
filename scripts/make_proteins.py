#!/usr/bin/env python
"""Translate PPanGGOLiN gene sequences into one protein FASTA per genome.

Proteins are written in genomic order (sorted by start in table/<genome>.tsv)
with PPanGGOLiN gene IDs, as MacSyFinder's ordered_replicon mode requires.
Translation uses NCBI table 11; the first codon is always written as M and the
trailing stop is removed. Internal stops (rare; fragmented genes) become X.
"""
import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def translate(nt):
    aa = str(Seq(nt).translate(table=11))
    aa = aa.rstrip("*")
    return ("M" + aa[1:]).replace("*", "X") if aa else aa


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", required=True, help="all_genes.fna from `ppanggolin fasta --genes all`")
    ap.add_argument("--tables", required=True, help="PPanGGOLiN table/ folder (one TSV per genome)")
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.genes, "fasta")}
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    for table in sorted(Path(args.tables).glob("*.tsv")):
        with open(table) as fh:
            genes = sorted(csv.DictReader(fh, delimiter="\t"),
                           key=lambda g: (g["contig"], int(g["start"])))
        n = 0
        with open(out / f"{table.stem}.faa", "w") as fh:
            for g in genes:
                if g["gene"] not in seqs:
                    continue  # non-CDS features (e.g. RNA genes)
                aa = translate(seqs[g["gene"]])
                fh.write(f">{g['gene']}\n")
                for i in range(0, len(aa), 60):
                    fh.write(aa[i:i + 60] + "\n")
                n += 1
        print(f"{table.stem}: {n} proteins")


if __name__ == "__main__":
    main()
