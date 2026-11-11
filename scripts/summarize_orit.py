#!/usr/bin/env python
"""Filter nhmmscan oriT hits and map them onto regions (milestone 2c).

A hit passes when:
  E-value <= --max-evalue (default 0.01), and
  aligned profile length >= --min-cov x profile length (default 0.7),
  where aligned profile length = hmm_to - hmm_from + 1.
The coverage filter removes the many short partial matches that oriT
profiles produce.

Outputs (in --outdir):
  orit_hits.tsv  every reported hit with coverage, pass flag and region_id
                 (empty region_id = hit lies outside the kept regions)
  orit.tsv       one row per kept region: orit (True/False), number of
                 passing hits, profiles hit (best E-value per profile)
"""
import argparse
from collections import defaultdict
from pathlib import Path

from common import find_region, read_tsv, region_segments, write_tsv

TBL_COLS = ["profile", "profile_acc", "contig", "contig_acc", "hmm_from", "hmm_to",
            "ali_from", "ali_to", "env_from", "env_to", "modlen", "strand",
            "evalue", "score", "bias"]


def profile_lengths(hmm_file):
    lengths, name = {}, None
    with open(hmm_file) as fh:
        for line in fh:
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("LENG"):
                lengths[name] = int(line.split()[1])
    return lengths


def read_tblout(path):
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=len(TBL_COLS))
            hits.append(dict(zip(TBL_COLS, fields[:len(TBL_COLS)])))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="regions/regions.tsv")
    ap.add_argument("--nhmmscan", required=True, help="folder with <genome>.tbl files")
    ap.add_argument("--hmm", required=True, help="oriT.hmm (for profile lengths)")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--min-cov", type=float, default=0.7,
                    help="minimum fraction of the oriT profile covered by the hit (default: 0.7)")
    ap.add_argument("--max-evalue", type=float, default=0.01,
                    help="maximum E-value (default: 0.01)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    regions = read_tsv(args.regions)
    segs = region_segments(regions)
    plen = profile_lengths(args.hmm)

    hits = []
    for tbl in sorted(Path(args.nhmmscan).glob("*.tbl")):
        for h in read_tblout(tbl):
            hmm_from, hmm_to = int(h["hmm_from"]), int(h["hmm_to"])
            a, b = int(h["ali_from"]), int(h["ali_to"])
            start, end = min(a, b), max(a, b)
            cov = (hmm_to - hmm_from + 1) / plen[h["profile"]]
            evalue = float(h["evalue"])
            hits.append({
                "genome": tbl.stem, "contig": h["contig"], "profile": h["profile"],
                "start": start, "end": end, "strand": h["strand"],
                "hmm_from": hmm_from, "hmm_to": hmm_to, "profile_len": plen[h["profile"]],
                "profile_cov": round(cov, 3), "evalue": evalue, "score": h["score"],
                "pass_evalue": evalue <= args.max_evalue,
                "pass_cov": cov >= args.min_cov,
                "region_id": find_region(segs, h["contig"], start, end),
            })
    for h in hits:
        h["pass"] = h["pass_evalue"] and h["pass_cov"]

    fields = ["genome", "contig", "profile", "start", "end", "strand", "hmm_from", "hmm_to",
              "profile_len", "profile_cov", "evalue", "score", "pass_evalue", "pass_cov",
              "pass", "region_id"]
    write_tsv(out / "orit_hits.tsv", sorted(hits, key=lambda h: (h["genome"], h["start"])), fields)

    per_region = defaultdict(list)
    for h in hits:
        if h["pass"] and h["region_id"]:
            per_region[h["region_id"]].append(h)
    rows = []
    for r in regions:
        rh = per_region.get(r["region_id"], [])
        best = {}
        for h in rh:
            best[h["profile"]] = min(best.get(h["profile"], 1.0), h["evalue"])
        rows.append({
            "region_id": r["region_id"],
            "orit": bool(rh),
            "n_orit_hits": len(rh),
            "orit_profiles": ";".join(f"{p}({e:.1e})" for p, e in sorted(best.items(), key=lambda x: x[1])),
            "orit_best_evalue": f"{min(best.values()):.1e}" if best else "",
        })
    write_tsv(out / "orit.tsv", rows)

    with open(out / "orit_filters.txt", "w") as fh:
        fh.write(f"max_evalue: {args.max_evalue}\nmin_profile_cov: {args.min_cov}\n")

    n_e = sum(h["pass_evalue"] for h in hits)
    n_pass = sum(h["pass"] for h in hits)
    n_in = sum(1 for h in hits if h["pass"] and h["region_id"])
    print(f"{len(hits)} hits reported; {n_e} pass E <= {args.max_evalue}; "
          f"{n_pass} also pass coverage >= {args.min_cov}; {n_in} of those in kept regions")
    print(f"regions with oriT: {sum(r['orit'] for r in rows)} of {len(rows)}")


if __name__ == "__main__":
    main()
