#!/usr/bin/env python
"""Search for direct repeats (DRs) at both ends of each region.

Integration of an ICE/IME by site-specific recombination duplicates the
attachment site, so the element ends up flanked by two identical copies
(attL and attR) in the same orientation.

For each region, two windows are compared:
  left  = [start - outside, start + inside)     around the region start
  right = (end - inside, end + outside]         around the region end
(--outside / --inside, default 3000 bp each; inside is capped at half the
region length; windows wrap around circular contigs.)

A DR candidate is an exact, same-strand match of >= --min-len bp (default 18)
present in both windows, extended to its maximal length. Candidates are kept if
  - they are not low-complexity (>= half of their 3-mers distinct), and
  - the DR is unique to its two flanking copies: the median genome-wide copy
    number (both strands) of its --min-len-mers is <= --max-copies (default 2).
    Using k-mers stops a REP element or other genome-wide repeat from passing
    just because a chance base extends it into a "unique" longer match.
The longest kept candidate is reported (ties: closest to the region ends).
Its outer coordinates give a refined element boundary: dr_left_start..dr_right_end.

Output (in --outdir): dr.tsv, one row per region; versions.txt
"""
import argparse
import csv
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import reverse_complement


def window(seq, start, end, circular):
    """1-based inclusive window [start, end].

    Returns (sequence, raw positions, contig positions). Raw positions are the
    unwrapped coordinates (may be < 1 or > len on circular contigs); contig
    positions are wrapped back onto the contig.
    """
    n = len(seq)
    raw = list(range(start, end + 1))
    if not circular:
        raw = [p for p in raw if 1 <= p <= n]
    pos = [((p - 1) % n) + 1 for p in raw]
    return "".join(seq[p - 1] for p in pos), raw, pos


def maximal_matches(a, b, k):
    """All maximal exact matches of length >= k between a and b: (len, i, j)."""
    index = {}
    for i in range(len(a) - k + 1):
        index.setdefault(a[i:i + k], []).append(i)
    found = []
    for j in range(len(b) - k + 1):
        for i in index.get(b[j:j + k], ()):
            if i > 0 and j > 0 and a[i - 1] == b[j - 1]:
                continue  # not left-maximal; reported from an earlier start
            length = k
            while i + length < len(a) and j + length < len(b) and a[i + length] == b[j + length]:
                length += 1
            found.append((length, i, j))
    return found


def low_complexity(s):
    trimers = {s[i:i + 3] for i in range(len(s) - 2)}
    return len(trimers) < 0.5 * (len(s) - 2)


def kmer_copies(seq, contigs, k):
    """Median genome-wide copy number of the k-mers of seq (both strands)."""
    counts = sorted(count_copies(seq[i:i + k], contigs) for i in range(len(seq) - k + 1))
    return counts[len(counts) // 2]


def count_copies(seq, contigs):
    """Occurrences (overlapping allowed) on both strands of all contigs."""
    pat = re.compile(f"(?=({re.escape(seq)}))")
    rc = reverse_complement(seq)
    pats = [pat] if rc == seq else [pat, re.compile(f"(?=({re.escape(rc)}))")]
    return sum(len(p.findall(c)) for c in contigs for p in pats)


def find_dr(contig_seq, all_contigs, start, end, args, circular):
    length = end - start + 1
    inside = min(args.inside, length // 2)
    left, lraw, lpos = window(contig_seq, start - args.outside, start + inside - 1, circular)
    right, rraw, rpos = window(contig_seq, end - inside + 1, end + args.outside, circular)

    seen, kept, n_raw = {}, [], 0
    for mlen, i, j in maximal_matches(left, right, args.min_len):
        n_raw += 1
        s = left[i:i + mlen]
        if low_complexity(s):
            continue
        if s not in seen:
            seen[s] = kmer_copies(s, all_contigs, args.min_len)
        if seen[s] > args.max_copies:
            continue
        # distance of each copy to its region end (for tie-breaking)
        dist = abs(lraw[i] - start) + abs(rraw[j + mlen - 1] - end)
        kept.append((mlen, -dist, i, j, s, seen[s]))
    if not kept:
        return {"dr": False, "n_dr_candidates": 0, "n_raw_matches": n_raw}

    mlen, _, i, j, s, copies = max(kept)
    ls, le = lpos[i], lpos[i + mlen - 1]
    rs, re_ = rpos[j], rpos[j + mlen - 1]
    return {
        "dr": True, "n_dr_candidates": len(kept), "n_raw_matches": n_raw,
        "dr_length": mlen, "dr_seq": s, "dr_genome_copies": copies,
        "dr_left_start": ls, "dr_left_end": le, "dr_right_start": rs, "dr_right_end": re_,
        # bp from the region end to the outer edge of each DR copy;
        # negative = outside the region (in the flank), positive = inside
        "dr_left_offset": lraw[i] - start,
        "dr_right_offset": end - rraw[j + mlen - 1],
        "element_start": ls, "element_end": re_,
    }


FIELDS = ["region_id", "dr", "dr_length", "dr_seq", "dr_genome_copies",
          "dr_left_start", "dr_left_end", "dr_right_start", "dr_right_end",
          "dr_left_offset", "dr_right_offset", "element_start", "element_end",
          "n_dr_candidates", "n_raw_matches"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="regions/regions.tsv")
    ap.add_argument("--fasta-list", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--outside", type=int, default=3000,
                    help="bp searched outside each region end (default: 3000)")
    ap.add_argument("--inside", type=int, default=3000,
                    help="bp searched inside each region end (default: 3000)")
    ap.add_argument("--min-len", type=int, default=18,
                    help="minimum DR length (default: 18; ~3.5%% background rate on the "
                         "test genomes with default windows, vs ~29%% at 15 bp)")
    ap.add_argument("--max-copies", type=int, default=2,
                    help="max genome-wide copies of a DR (default: 2)")
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
        contigs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(genomes[genome], "fasta")}
        all_seqs = list(contigs.values())
        for r in (r for r in regions if r["genome"] == genome):
            segs = [tuple(map(int, p.split(".."))) for p in r["coordinates"].split(",")]
            start, end = segs[0][0], segs[-1][1]
            if end < start:  # region wraps the origin: express end past the contig length
                end += len(contigs[r["contig"]])
            res = find_dr(contigs[r["contig"]], all_seqs, start, end, args,
                          r["contig"] in circular)
            rows.append({"region_id": r["region_id"], **res})

    with open(out / "dr.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    with open(out / "versions.txt", "w") as fh:
        fh.write(f"date: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n"
                 f"python: {sys.version.split()[0]}\n"
                 f"outside_bp: {args.outside}\ninside_bp: {args.inside}\n"
                 f"min_len: {args.min_len}\nmax_copies: {args.max_copies}\n"
                 f"seed: not applicable (deterministic)\n")
    print(f"DRs found in {sum(r['dr'] for r in rows)} of {len(rows)} regions -> {out / 'dr.tsv'}")


if __name__ == "__main__":
    main()
