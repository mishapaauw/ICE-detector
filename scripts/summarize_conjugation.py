#!/usr/bin/env python
"""Map MacSyFinder/CONJScan hits onto PPanGGOLiN regions (milestone 2a + 2b).

Inputs:
  --regions       regions/regions.tsv (regions kept after merging and --min-orfs)
  --region-genes  regions/region_genes.tsv: gene -> region (incl. merge-gap genes)
  --tables        PPanGGOLiN table/ folder: gene -> RGP membership for all RGPs
  --macsyfinder   folder with one MacSyFinder output folder per genome
Outputs (in --outdir):
  conjugation_hits.tsv  every hit of the best solution, with the kept region
                        (region_id) and the PPanGGOLiN RGP (rgp_id, any size) it lies in
  conjugation.tsv       one row per kept region:
    conj_system      True if a T4SS system has >= --min-frac of its hits in the region
    conj_systems     e.g. "T4SS_typeG:barcode71_rotated_T4SS_typeG_6(21/21)"
                     = system type : system id (hits in region / hits in system)
    dconj_system     same for decayed conjugation systems (dCONJ_*)
    relaxase         True if any relaxase (MOB*) hit lies in the region
    relaxase_types   relaxase families found in the region, e.g. "MOBH,MOBQ"
    t4cp, virb4      True if a coupling protein / VirB4 ATPase hit lies in the region
"""
import argparse
from collections import defaultdict
from pathlib import Path

from common import read_tsv, write_tsv


def gene_to_rgp(tables_dir):
    mapping = {}
    for table in Path(tables_dir).glob("*.tsv"):
        for g in read_tsv(table):
            if g["RGP"]:
                mapping[g["gene"]] = g["RGP"]
    return mapping


def load_hits(msf_dir, region_of, rgp_of):
    hits = []
    for best in sorted(Path(msf_dir).glob("*/best_solution.tsv")):
        for h in read_tsv(best, comment="#"):
            if not h.get("hit_id"):
                continue
            model = h["model_fqn"].split("/")[-1]  # e.g. T4SS_typeG, dCONJ_typeF, MOB
            hits.append({
                "genome": best.parent.name, "hit_id": h["hit_id"],
                "gene_name": h["gene_name"], "model": model, "sys_id": h["sys_id"],
                "hit_status": h["hit_status"], "hit_i_eval": h["hit_i_eval"],
                "region_id": region_of.get(h["hit_id"], ""),
                "rgp_id": rgp_of.get(h["hit_id"], ""),
            })
    return hits


def relaxase_family(gene_name):
    # T4SS_MOBQ -> MOBQ ; exchangeable names follow the same pattern
    return gene_name.split("_")[-1] if "_MOB" in gene_name else None


def summarize(regions, hits, min_frac):
    sys_size = defaultdict(int)
    for h in hits:
        sys_size[h["sys_id"]] += 1

    per_region = defaultdict(list)
    for h in hits:
        if h["region_id"]:
            per_region[h["region_id"]].append(h)

    rows = []
    for reg in regions:
        rid = reg["region_id"]
        rh = per_region.get(rid, [])
        in_region = defaultdict(int)
        model_of = {}
        for h in rh:
            in_region[h["sys_id"]] += 1
            model_of[h["sys_id"]] = h["model"]

        def systems(prefix):
            out = []
            for sid, n in sorted(in_region.items()):
                if model_of[sid].startswith(prefix):
                    out.append((model_of[sid], sid, n, sys_size[sid]))
            return out

        t4ss = systems("T4SS")
        dconj = systems("dCONJ")
        fmt = lambda s: ";".join(f"{m}:{sid}({n}/{tot})" for m, sid, n, tot in s)
        mobs = sorted({f for h in rh if (f := relaxase_family(h["gene_name"]))})
        rows.append({
            "region_id": rid,
            "conj_system": any(n / tot >= min_frac for _, _, n, tot in t4ss),
            "conj_systems": fmt(t4ss),
            "dconj_system": any(n / tot >= min_frac for _, _, n, tot in dconj),
            "dconj_systems": fmt(dconj),
            "relaxase": bool(mobs),
            "relaxase_types": ",".join(mobs),
            "t4cp": any("t4cp" in h["gene_name"] or "tcpA" in h["gene_name"] for h in rh),
            "virb4": any(h["gene_name"] in ("T4SS_virb4", "T4SS_I_traU") for h in rh),
            "n_conj_hits": len(rh),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--region-genes", required=True)
    ap.add_argument("--tables", required=True)
    ap.add_argument("--macsyfinder", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="fraction of a system's hits that must lie in a region "
                         "to call conj_system=True (default: 0.5)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    regions = read_tsv(args.regions)
    region_of = {g["gene"]: g["region_id"] for g in read_tsv(args.region_genes)}
    hits = load_hits(args.macsyfinder, region_of, gene_to_rgp(args.tables))
    rows = summarize(regions, hits, args.min_frac)

    write_tsv(out / "conjugation_hits.tsv", hits, list(hits[0]) if hits else ["hit_id"])
    write_tsv(out / "conjugation.tsv", rows, list(rows[0]))

    kept = {r["region_id"] for r in regions}
    n_sys = len({h["sys_id"] for h in hits})
    outside = {h["sys_id"] for h in hits if h["region_id"] not in kept}
    print(f"{len(hits)} hits in {n_sys} systems; {len(outside)} systems have hits "
          f"outside the kept regions (see conjugation_hits.tsv)")
    print(f"regions with conj_system: {sum(r['conj_system'] for r in rows)}, "
          f"relaxase: {sum(r['relaxase'] for r in rows)} (of {len(rows)})")


if __name__ == "__main__":
    main()
