#!/usr/bin/env python
"""Predict tRNA and tmRNA genes on every chromosome with ARAGORN (pyaragorn).

Used twice downstream: extract_regions.py does not merge RGPs across a tRNA,
and trna_proximity.py flags regions with a tRNA at one end.

Output (in --outdir):
  trnas.tsv     genome, contig, contig_length, circular, type, name
                (e.g. tRNA-Leu(cag), tmRNA), start, end (1-based, start <= end), strand
  versions.txt
"""
import argparse
from pathlib import Path

import pyaragorn
from Bio import SeqIO

from common import read_fasta_list, write_tsv, write_versions

FIELDS = ["genome", "contig", "contig_length", "circular", "type", "name", "start", "end",
          "strand"]


def predict(seq, circular):
    """tRNA/tmRNA genes of one contig: list of dicts with 1-based start <= end."""
    finder = pyaragorn.RNAFinder(translation_table=11, trna=True, tmrna=True, linear=not circular)
    genes = []
    for g in finder.find_rna(seq):
        start, end = sorted((g.begin, g.end))
        if g.type == "tRNA":
            name = f"tRNA-{g.amino_acid}({g.anticodon.lower()})"
        else:
            name = g.type
        genes.append({"type": g.type, "name": name, "start": start, "end": end,
                      "strand": "+" if g.strand > 0 else "-"})
    return genes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta-list", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    genomes, circular = read_fasta_list(args.fasta_list)
    rows = []
    for genome in sorted(genomes):
        for rec in SeqIO.parse(genomes[genome], "fasta"):
            seq, circ = str(rec.seq).upper(), rec.id in circular
            rows += [{"genome": genome, "contig": rec.id, "contig_length": len(seq),
                      "circular": circ, **g} for g in predict(seq, circ)]

    write_tsv(out / "trnas.tsv", rows, FIELDS)
    write_versions(out / "versions.txt", pyaragorn=pyaragorn.__version__,
                   aragorn_settings="genetic code 11, tRNA + tmRNA",
                   seed="not applicable (deterministic)")
    print(f"{len(rows)} tRNA/tmRNA genes in {len(genomes)} genomes -> {out / 'trnas.tsv'}")


if __name__ == "__main__":
    main()
