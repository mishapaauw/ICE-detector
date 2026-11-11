#!/usr/bin/env python
"""tRNA / tmRNA genes next to region ends (integration hotspots).

Many ICEs and other genomic islands integrate into the 3' end of a tRNA or tmRNA
gene, so a tRNA at one end of a region is supporting evidence for an integrated
element. tRNAs/tmRNAs come from predict_trnas.py (ARAGORN via pyaragorn).

Per region:
  trna_adjacent   True if a tRNA/tmRNA lies within --dist bp of a region end
                  (inside or outside the region), or overlaps a direct-repeat copy
  trna_name       the nearest such gene, e.g. "tRNA-Leu(cag)" or "tmRNA"
  trna_end        region end it is near: left (start) or right (end)
  trna_dist       bp between the gene and that end (0 = overlaps the end)
  trna_dr_overlap True if the gene overlaps a DR copy (att site inside the tRNA)

Inputs: regions.tsv, trnas.tsv (predict_trnas.py), optional dr/dr.tsv.
Outputs (in --outdir): trna.tsv (per region); settings appended to versions.txt
"""
import argparse
from pathlib import Path

from common import parse_coordinates, read_tsv, write_tsv


def point_distance(p, start, end, n, circular):
    """bp between position p and the interval start..end (0 if inside)."""
    if start <= p <= end:
        return 0
    d = min(abs(p - start), abs(p - end))
    return min(d, n - d) if circular else d


def overlaps(a1, a2, b1, b2):
    return a1 <= b2 and b1 <= a2


def nearest_trna(genes, left, right, n, circular, max_dist, dr=None):
    """Closest tRNA to either region end within max_dist, or one overlapping a DR copy."""
    dr_copies = []
    if dr and str(dr.get("dr")) == "True":
        dr_copies = [(int(dr["dr_left_start"]), int(dr["dr_left_end"])),
                     (int(dr["dr_right_start"]), int(dr["dr_right_end"]))]
    best = None
    for g in genes:
        on_dr = any(overlaps(g["start"], g["end"], s, e) for s, e in dr_copies)
        for side, p in (("left", left), ("right", right)):
            d = point_distance(p, g["start"], g["end"], n, circular)
            if d <= max_dist or on_dr:
                cand = (d, side, g, on_dr)
                if best is None or (cand[3], -cand[0]) > (best[3], -best[0]):
                    best = cand  # prefer DR overlap, then the shortest distance
    if best is None:
        return {"trna_adjacent": False, "trna_name": "", "trna_end": "", "trna_dist": "",
                "trna_dr_overlap": False}
    d, side, g, on_dr = best
    return {"trna_adjacent": True, "trna_name": g["name"], "trna_end": side, "trna_dist": d,
            "trna_dr_overlap": on_dr}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="regions/regions.tsv")
    ap.add_argument("--trnas", required=True, help="trnas.tsv from predict_trnas.py")
    ap.add_argument("--dr", help="dr/dr.tsv (optional): also count tRNAs overlapping a DR copy")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--dist", type=int, default=1000,
                    help="max bp between a tRNA/tmRNA and a region end (default: 1000)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    regions = read_tsv(args.regions)
    dr = {r["region_id"]: r for r in read_tsv(args.dr)} if args.dr else {}
    genes, contig_info = {}, {}
    for g in read_tsv(args.trnas):
        g["start"], g["end"] = int(g["start"]), int(g["end"])
        genes.setdefault(g["contig"], []).append(g)
        contig_info[g["contig"]] = (int(g["contig_length"]), g["circular"] == "True")

    rows = []
    for r in regions:
        segs = parse_coordinates(r["coordinates"])
        c = r["contig"]
        if c in contig_info:
            n, circ = contig_info[c]
            res = nearest_trna(genes[c], segs[0][0], segs[-1][1], n, circ, args.dist,
                               dr.get(r["region_id"]))
        else:  # no tRNA predicted on this contig
            res = nearest_trna([], 0, 0, 1, False, args.dist)
        rows.append({"region_id": r["region_id"], **res})

    write_tsv(out / "trna.tsv", rows, ["region_id", "trna_adjacent", "trna_name", "trna_end",
                                       "trna_dist", "trna_dr_overlap"])
    with open(out / "versions.txt", "a") as fh:
        fh.write(f"max_dist_bp: {args.dist}\ndr_table: {args.dr or 'not used'}\n")
    print(f"regions with an adjacent tRNA: "
          f"{sum(r['trna_adjacent'] for r in rows)} of {len(rows)} -> {out / 'trna.tsv'}")


if __name__ == "__main__":
    main()
