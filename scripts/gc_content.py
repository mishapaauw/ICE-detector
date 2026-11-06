#!/usr/bin/env python
"""GC content of each region compared to its genome (milestone 2d).

For each region:
  region_gc    GC fraction of the region (over A/C/G/T only)
  genome_gc    GC fraction of the whole genome
  gc_diff      region_gc - genome_gc
  gc_zscore    z-score of region_gc against all same-length windows of the
               same contig (step --step bp; windows wrap around circular contigs).
               Strongly negative/positive = atypical composition, as expected
               for recently acquired elements.
  flank_up_gc, flank_down_gc
               GC fraction of the --flank bp (default 35 kb) directly upstream
               (before start) and downstream (after end) of the region. Flanks
               wrap around circular contigs and are truncated at linear contig ends.
  flank_gc     GC of both flanks together
  gc_diff_flank  region_gc - flank_gc (local contrast, independent of
               genome-wide GC variation)
"""
import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio import SeqIO


def counts(seq):
    arr = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    gc = (arr == ord("G")) | (arr == ord("C"))
    acgt = gc | (arr == ord("A")) | (arr == ord("T"))
    return gc.astype(np.int64), acgt.astype(np.int64)


def window_gc(gc, acgt, length, step, circular):
    """GC fraction of all windows of `length` bp, every `step` bp."""
    n = len(gc)
    if circular:
        gc = np.concatenate([gc, gc[:length]])
        acgt = np.concatenate([acgt, acgt[:length]])
        starts = np.arange(0, n, step)
    else:
        starts = np.arange(0, max(n - length, 0) + 1, step)
    cgc = np.concatenate([[0], np.cumsum(gc)])
    cat = np.concatenate([[0], np.cumsum(acgt)])
    g = cgc[starts + length] - cgc[starts]
    t = cat[starts + length] - cat[starts]
    ok = t > 0
    return g[ok] / t[ok]


def flank_gc(gc, acgt, start, end, flank, circular):
    """GC of `flank` bp before `start` and after `end` (1-based, inclusive).

    Returns (up_gc, down_gc, both_gc); None where a flank has no bases.
    """
    n = len(gc)
    up = np.arange(start - 1 - flank, start - 1)
    down = np.arange(end, end + flank)
    if circular:
        up, down = up % n, down % n
    else:
        up, down = up[up >= 0], down[down < n]

    def frac(idx):
        t = int(acgt[idx].sum())
        return int(gc[idx].sum()) / t if t else None

    both = np.concatenate([up, down])
    return frac(up), frac(down), frac(both)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="regions/regions.tsv")
    ap.add_argument("--fasta-list", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--step", type=int, default=1000, help="window step in bp (default: 1000)")
    ap.add_argument("--flank", type=int, default=35000,
                    help="flank size up- and downstream of each region in bp (default: 35000)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.regions) as fh:
        regions = list(csv.DictReader(fh, delimiter="\t"))
    genomes, circular = {}, set()
    with open(args.fasta_list) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row:
                genomes[row[0]] = row[1]
                circular.update(row[2:])

    rows = []
    for genome in sorted({r["genome"] for r in regions}):
        contigs = {rec.id: str(rec.seq) for rec in SeqIO.parse(genomes[genome], "fasta")}
        cnt = {c: counts(s) for c, s in contigs.items()}
        g_gc = sum(int(g.sum()) for g, _ in cnt.values())
        g_at = sum(int(a.sum()) for _, a in cnt.values())
        genome_gc = g_gc / g_at
        for r in (r for r in regions if r["genome"] == genome):
            gc, acgt = cnt[r["contig"]]
            rg = ra = 0
            for part in r["coordinates"].split(","):
                s, e = map(int, part.split(".."))
                rg += int(gc[s - 1:e].sum())
                ra += int(acgt[s - 1:e].sum())
            region_gc = rg / ra
            segs = [tuple(map(int, p.split(".."))) for p in r["coordinates"].split(",")]
            up, down, both = flank_gc(gc, acgt, segs[0][0], segs[-1][1], args.flank,
                                      r["contig"] in circular)
            rnd = lambda x: round(x, 4) if x is not None else ""
            win = window_gc(gc, acgt, int(r["length"]), args.step, r["contig"] in circular)
            sd = win.std()
            rows.append({
                "region_id": r["region_id"],
                "region_gc": round(region_gc, 4),
                "genome_gc": round(genome_gc, 4),
                "gc_diff": round(region_gc - genome_gc, 4),
                "gc_zscore": round((region_gc - win.mean()) / sd, 2) if sd > 0 else "",
                "flank_up_gc": rnd(up),
                "flank_down_gc": rnd(down),
                "flank_gc": rnd(both),
                "gc_diff_flank": rnd(region_gc - both) if both is not None else "",
            })

    with open(out / "gc.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(out / "versions.txt", "w") as fh:
        fh.write(f"date: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n"
                 f"python: {sys.version.split()[0]}\nnumpy: {np.__version__}\n"
                 f"window_step: {args.step}\nflank_bp: {args.flank}\nseed: not applicable (deterministic)\n")
    print(f"GC computed for {len(rows)} regions -> {out / 'gc.tsv'}")


if __name__ == "__main__":
    main()
