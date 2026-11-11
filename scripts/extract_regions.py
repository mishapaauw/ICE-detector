#!/usr/bin/env python
"""Extract PPanGGOLiN regions of genomic plasticity (RGPs) as FASTA + table.

Neighbouring RGPs on the same contig that are at most --merge-gap genes apart are
merged into one region first: one element is often split into several RGPs by a
gene that is persistent in the pangenome (e.g. a cargo gene shared with the core).
The --min-orfs filter is applied to the merged region, so the pieces are not lost.
RGPs are never merged across a tRNA/tmRNA gene (with --trnas): a tRNA between two
RGPs usually marks the integration site of a separate element. RGPs that wrap the
origin of a circular contig are never merged.

Inputs (from a PPanGGOLiN output folder):
  regions_of_genomic_plasticity.tsv, spots.tsv, table/<genome>.tsv
Outputs (in --outdir):
  regions.tsv        one row per region with >= --min-orfs genes; region_id = ID of
                     its largest RGP, rgp_ids = all RGPs merged into it (genomic order)
  region_genes.tsv   genes inside those regions, including genes in merge gaps
  fasta/<region_id>.fasta  genomic sequence of each region, cut from the input genomes
  versions.txt
"""
import argparse
import sys
from pathlib import Path

import Bio
from Bio import SeqIO

from common import parse_coordinates, read_fasta_list, read_tsv, write_tsv, write_versions

REGION_COLS = ["region_id", "spot_id", "genome", "contig", "start", "end",
               "length", "n_orfs", "coordinates", "contig_border", "whole_contig",
               "rgp_ids", "n_merged"]
GENE_COLS = ["region_id", "gene", "contig", "start", "stop", "strand", "family", "partition"]


def cut_sequence(seq, segments):
    """Concatenate 1-based inclusive segments of seq, in the given order."""
    return "".join(str(seq[s - 1:e]) for s, e in segments)


def gene_positions(genes):
    """gene -> (contig, index in genomic order on that contig)."""
    pos, by_contig = {}, {}
    for g in genes:
        by_contig.setdefault(g["contig"], []).append(g)
    for contig, gs in by_contig.items():
        for i, g in enumerate(sorted(gs, key=lambda g: int(g["start"]))):
            pos[g["gene"]] = (contig, i)
    return pos


def trna_between(trnas, contig, left_end, right_start):
    """True if a tRNA on contig overlaps the gap left_end < x < right_start."""
    return any(s < right_start and e > left_end for s, e in trnas.get(contig, ()))


def merge_rgps(rgps, pos, gap, trnas=None):
    """Group RGPs of one genome that are <= gap genes apart on the same contig,
    with no tRNA/tmRNA in between.

    rgps: RGP rows (first_gene, last_gene, start, stop, coordinates, ...);
    pos: gene_positions(); trnas: {contig: [(start, end), ...]} or None.
    Returns lists of RGPs, each in genomic order, ordered by the input position
    of their first member (so gap=0 keeps the input order).
    """
    order = {r["region"]: i for i, r in enumerate(rgps)}
    linear = [r for r in rgps if "," not in r["coordinates"]]  # wrapping RGPs stay alone
    groups = [[r] for r in rgps if "," in r["coordinates"]]
    if gap > 0:
        linear.sort(key=lambda r: pos[r["first_gene"]])
        for r in linear:
            prev = groups[-1][-1] if groups and "," not in groups[-1][-1]["coordinates"] else None
            if prev is not None:
                (c1, i1), (c2, i2) = pos[prev["last_gene"]], pos[r["first_gene"]]
                if (c1 == c2 and 0 <= i2 - i1 - 1 <= gap and
                        not (trnas and trna_between(trnas, c1, int(prev["stop"]),
                                                    int(r["start"])))):
                    groups[-1].append(r)
                    continue
            groups.append([r])
    else:
        groups += [[r] for r in linear]
    return sorted(groups, key=lambda g: order[g[0]["region"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppanggolin", required=True, help="PPanGGOLiN output folder")
    ap.add_argument("--fasta-list", required=True, help="fasta list given to PPanGGOLiN")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--min-orfs", type=int, default=10,
                    help="keep regions with at least this many genes (default: 10)")
    ap.add_argument("--trnas", help="trnas.tsv (predict_trnas.py): do not merge across a tRNA")
    ap.add_argument("--merge-gap", type=int, default=3,
                    help="merge neighbouring RGPs at most this many genes apart "
                         "(default: 3; 0 = no merging)")
    args = ap.parse_args()

    ppg = Path(args.ppanggolin)
    out = Path(args.outdir)
    (out / "fasta").mkdir(parents=True, exist_ok=True)

    spot_of = {r["rgp_id"]: r["spot_id"] for r in read_tsv(ppg / "spots.tsv")}
    all_rgps = read_tsv(ppg / "regions_of_genomic_plasticity.tsv")
    genome_paths, _ = read_fasta_list(args.fasta_list)
    trnas = {}
    for t in (read_tsv(args.trnas) if args.trnas else []):
        trnas.setdefault(t["genome"], {}).setdefault(t["contig"], []).append(
            (int(t["start"]), int(t["end"])))
    by_genome = {}
    for r in all_rgps:
        by_genome.setdefault(r["genome"], []).append(r)

    rows, gene_rows, n_merged_regions = [], [], 0
    for genome, rgps in sorted(by_genome.items()):
        genes = read_tsv(ppg / "table" / f"{genome}.tsv")
        pos = gene_positions(genes)
        contigs = None
        for group in merge_rgps(rgps, pos, args.merge_gap, trnas.get(genome)):
            first, last = group[0], group[-1]
            # merged regions keep the ID of their largest RGP (stable between settings)
            rid = max(group, key=lambda r: int(r["genes"]))["region"]
            contig = first["contig"]
            if len(group) == 1:
                members = {rid}
                region_genes = [g for g in genes if g["RGP"] == rid]
                n_orfs, coords = int(first["genes"]), first["coordinates"]
                start, end, length = first["start"], first["stop"], int(first["length"])
            else:
                lo, hi = pos[first["first_gene"]][1], pos[last["last_gene"]][1]
                region_genes = [g for g in genes if pos[g["gene"]][0] == contig
                                and lo <= pos[g["gene"]][1] <= hi]
                n_orfs = len(region_genes)
                start, end = first["start"], last["stop"]
                coords, length = f"{start}..{end}", int(end) - int(start) + 1
            if n_orfs < args.min_orfs:
                continue
            n_merged_regions += len(group) > 1
            if contigs is None:
                contigs = {rec.id: rec.seq for rec in SeqIO.parse(genome_paths[genome], "fasta")}
            seq = cut_sequence(contigs[contig], parse_coordinates(coords))
            if len(seq) != length:
                sys.exit(f"{rid}: extracted {len(seq)} bp but expected {length}")
            spots = list(dict.fromkeys(spot_of[r["region"]] for r in group
                                       if r["region"] in spot_of)) or ["no_spot"]
            spot = ";".join(spots)
            with open(out / "fasta" / f"{rid}.fasta", "w") as fh:
                fh.write(f">{rid} genome={genome} contig={contig} "
                         f"coords={coords} spot={spot}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
            rows.append({
                "region_id": rid, "spot_id": spot, "genome": genome,
                "contig": contig, "start": start, "end": end,
                "length": length, "n_orfs": n_orfs, "coordinates": coords,
                "contig_border": any(r["contigBorder"] == "True" for r in group),
                "whole_contig": any(r["wholeContig"] == "True" for r in group),
                "rgp_ids": ";".join(r["region"] for r in group), "n_merged": len(group),
            })
            gene_rows += [{"region_id": rid, **{k: g[k] for k in GENE_COLS[1:]}}
                          for g in region_genes]
    if not rows:
        sys.exit(f"No regions with >= {args.min_orfs} genes.")

    write_tsv(out / "regions.tsv", rows, REGION_COLS)
    write_tsv(out / "region_genes.tsv", gene_rows, GENE_COLS)
    write_versions(out / "versions.txt", biopython=Bio.__version__,
                   min_orfs=args.min_orfs, merge_gap=args.merge_gap,
                   merge_blocked_by_trna=bool(args.trnas),
                   ppanggolin_dir=ppg.resolve())

    print(f"{len(all_rgps)} RGPs total; {len(rows)} regions with >= {args.min_orfs} ORFs, "
          f"{n_merged_regions} of them merged from neighbouring RGPs "
          f"(--merge-gap {args.merge_gap}) -> {out}/regions.tsv, {out}/fasta/")


if __name__ == "__main__":
    main()
