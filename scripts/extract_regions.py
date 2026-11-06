#!/usr/bin/env python
"""Extract PPanGGOLiN regions of genomic plasticity (RGPs) as FASTA + table.

Inputs (from a PPanGGOLiN output folder):
  regions_of_genomic_plasticity.tsv, spots.tsv, table/<genome>.tsv
Outputs (in --outdir):
  regions.tsv        one row per RGP with >= --min-orfs genes
  region_genes.tsv   genes inside those RGPs (coordinates, family, partition)
  fasta/<rgp>.fasta  genomic sequence of each RGP, cut from the input genomes
  versions.txt
"""
import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import SeqIO

REGION_COLS = ["region_id", "spot_id", "genome", "contig", "start", "end",
               "length", "n_orfs", "coordinates", "contig_border", "whole_contig"]


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_coordinates(coords):
    """'100..200' or '5000..5100,1..50' (wraps origin) -> [(100, 200), ...]."""
    segments = []
    for part in coords.split(","):
        start, end = part.split("..")
        segments.append((int(start), int(end)))
    return segments


def cut_sequence(seq, segments):
    """Concatenate 1-based inclusive segments of seq, in the given order."""
    return "".join(str(seq[s - 1:e]) for s, e in segments)


def load_genome_paths(fasta_list):
    with open(fasta_list) as handle:
        return {row[0]: row[1] for row in csv.reader(handle, delimiter="\t") if row}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppanggolin", required=True, help="PPanGGOLiN output folder")
    ap.add_argument("--fasta-list", required=True, help="fasta list given to PPanGGOLiN")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--min-orfs", type=int, default=10,
                    help="keep RGPs with at least this many genes (default: 10)")
    args = ap.parse_args()

    ppg = Path(args.ppanggolin)
    out = Path(args.outdir)
    (out / "fasta").mkdir(parents=True, exist_ok=True)

    spot_of = {r["rgp_id"]: r["spot_id"] for r in read_tsv(ppg / "spots.tsv")}
    all_rgps = read_tsv(ppg / "regions_of_genomic_plasticity.tsv")
    kept = [r for r in all_rgps if int(r["genes"]) >= args.min_orfs]
    if not kept:
        sys.exit(f"No RGPs with >= {args.min_orfs} genes.")

    genome_paths = load_genome_paths(args.fasta_list)
    by_genome = {}
    for r in kept:
        by_genome.setdefault(r["genome"], []).append(r)

    rows, gene_rows = [], []
    for genome, rgps in sorted(by_genome.items()):
        contigs = {rec.id: rec.seq for rec in SeqIO.parse(genome_paths[genome], "fasta")}
        genes = read_tsv(ppg / "table" / f"{genome}.tsv")
        for r in rgps:
            rid = r["region"]
            seq = cut_sequence(contigs[r["contig"]], parse_coordinates(r["coordinates"]))
            if len(seq) != int(r["length"]):
                sys.exit(f"{rid}: extracted {len(seq)} bp but PPanGGOLiN says {r['length']}")
            spot = spot_of.get(rid, "no_spot")
            with open(out / "fasta" / f"{rid}.fasta", "w") as fh:
                fh.write(f">{rid} genome={genome} contig={r['contig']} "
                         f"coords={r['coordinates']} spot={spot}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
            rows.append({
                "region_id": rid, "spot_id": spot, "genome": genome,
                "contig": r["contig"], "start": r["start"], "end": r["stop"],
                "length": r["length"], "n_orfs": r["genes"],
                "coordinates": r["coordinates"],
                "contig_border": r["contigBorder"], "whole_contig": r["wholeContig"],
            })
            gene_rows += [{"region_id": rid, **{k: g[k] for k in
                           ("gene", "contig", "start", "stop", "strand", "family", "partition")}}
                          for g in genes if g["RGP"] == rid]

    with open(out / "regions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REGION_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(out / "region_genes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gene_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)
    with open(out / "versions.txt", "w") as fh:
        fh.write(f"date: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n"
                 f"python: {sys.version.split()[0]}\nbiopython: {Bio.__version__}\n"
                 f"min_orfs: {args.min_orfs}\nppanggolin_dir: {ppg.resolve()}\n")

    print(f"{len(all_rgps)} RGPs total, {len(rows)} with >= {args.min_orfs} ORFs "
          f"-> {out}/regions.tsv, {out}/fasta/")


if __name__ == "__main__":
    main()
