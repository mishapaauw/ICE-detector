#!/usr/bin/env python
"""Map integrase / recombinase hits (hmmsearch --cut_ga) onto regions.

Inputs:
  --regions       regions/regions.tsv
  --region-genes  regions/region_genes.tsv (genes in each region, with coordinates)
  --hmmsearch     folder with one <genome>.domtbl per genome (detect_integrase.sh)
Outputs (in --outdir):
  integrase_hits.tsv  one row per protein x profile (best domain), with its region
  integrase.tsv       one row per region:
    integrase          True if a tyrosine (Tyr) or serine (Ser) recombinase lies in the region
    integrase_types    e.g. "Tyr" or "Ser,Tyr"
    n_integrases       number of Tyr/Ser recombinase genes in the region
    integrase_end_dist bp from the closest Tyr/Ser recombinase gene to the nearest region end
                       (integrases usually sit next to one attachment site)
    dde_transposase    True if a DDE transposase (Pfam rve) lies in the region; reported
                       only, not used as evidence: IS elements carry it everywhere
"""
import argparse
from collections import defaultdict
from pathlib import Path

from common import parse_coordinates, read_tsv, write_tsv

# Pfam accession (without version) -> role
ROLES = {"PF00589": "Tyr", "PF00239": "Ser", "PF07508": "Ser", "PF00665": "DDE"}
INTEGRASE_ROLES = {"Tyr", "Ser"}


def read_domtbl(path):
    """Best domain (lowest i-Evalue) per (protein, profile) from a --domtblout file."""
    best = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            hit = {"gene": f[0], "profile": f[3], "accession": f[4].split(".")[0],
                   "evalue": float(f[6]), "score": float(f[7]), "dom_ievalue": float(f[12]),
                   "hmm_from": int(f[15]), "hmm_to": int(f[16]), "profile_len": int(f[5])}
            key = (hit["gene"], hit["profile"])
            if key not in best or hit["dom_ievalue"] < best[key]["dom_ievalue"]:
                best[key] = hit
    return list(best.values())


def offset_in_region(pos, segs):
    """0-based offset of contig position `pos` from the region start (segments in order)."""
    done = 0
    for s, e in segs:
        if s <= pos <= e:
            return done + pos - s
        done += e - s + 1
    return None


def end_distance(gene, segs, length):
    """bp from a gene to the nearest region end (0 = touches an end)."""
    a = offset_in_region(int(gene["start"]), segs)
    b = offset_in_region(int(gene["stop"]), segs)
    if a is None or b is None:  # gene sticks out of the region
        return 0
    return max(0, min(min(a, b), length - 1 - max(a, b)))


def summarize(regions, region_genes, hits):
    """Add role/region to hits (in place) and return one summary row per region."""
    gene_info = {g["gene"]: g for g in region_genes}
    for h in hits:
        h["role"] = ROLES.get(h["accession"], "other")
        h["region_id"] = gene_info[h["gene"]]["region_id"] if h["gene"] in gene_info else ""

    per_region = defaultdict(list)
    for h in hits:
        if h["region_id"]:
            per_region[h["region_id"]].append(h)

    rows = []
    for r in regions:
        rh = per_region.get(r["region_id"], [])
        integ = [h for h in rh if h["role"] in INTEGRASE_ROLES]
        genes = {h["gene"] for h in integ}
        segs, length = parse_coordinates(r["coordinates"]), int(r["length"])
        dists = [end_distance(gene_info[g], segs, length) for g in genes]
        rows.append({
            "region_id": r["region_id"],
            "integrase": bool(integ),
            "integrase_types": ",".join(sorted({h["role"] for h in integ})),
            "n_integrases": len(genes),
            "integrase_end_dist": min(dists) if dists else "",
            "dde_transposase": any(h["role"] == "DDE" for h in rh),
        })
    return rows


HIT_FIELDS = ["genome", "gene", "profile", "accession", "role", "evalue", "score",
              "dom_ievalue", "hmm_from", "hmm_to", "profile_len", "region_id"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--region-genes", required=True)
    ap.add_argument("--hmmsearch", required=True, help="folder with <genome>.domtbl files")
    ap.add_argument("-o", "--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    hits = []
    for tbl in sorted(Path(args.hmmsearch).glob("*.domtbl")):
        for h in read_domtbl(tbl):
            hits.append({"genome": tbl.stem, **h})
    regions = read_tsv(args.regions)
    rows = summarize(regions, read_tsv(args.region_genes), hits)

    write_tsv(out / "integrase_hits.tsv", sorted(hits, key=lambda h: (h["genome"], h["gene"])),
              HIT_FIELDS)
    write_tsv(out / "integrase.tsv", rows)

    n = defaultdict(int)
    for h in hits:
        n[h["role"]] += 1
    print("hits: " + ", ".join(f"{k}={v}" for k, v in sorted(n.items())))
    print(f"regions with integrase (Tyr/Ser): {sum(r['integrase'] for r in rows)}, "
          f"with DDE transposase: {sum(r['dde_transposase'] for r in rows)} (of {len(rows)})")


if __name__ == "__main__":
    main()
