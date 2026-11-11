#!/usr/bin/env python
"""Search for direct repeats (DRs) at both ends of each region.

Integration of an ICE/IME by site-specific recombination duplicates the
attachment site, so the element ends up flanked by two identical copies
(attL and attR) in the same orientation.

For each region, two windows are compared:
  left  = [start - outside, start + inside)     around the region start
  right = (end - inside, end + outside]         around the region end
(--outside / --inside, default 5000 bp each; inside is capped at half the
region length; windows wrap around circular contigs.)

A DR candidate is an exact, same-strand match of >= --min-len bp (default 15)
present in both windows, extended to its maximal length. Candidates are kept if
  - they are not low-complexity (>= half of their 3-mers distinct), and
  - the DR is unique to its two flanking copies: it occurs at most --max-copies
    times in the genome (both strands; default 2), and versions trimmed by up to
    2 bp at either end occur at most --max-copies + 1 times. The trimmed check
    stops a REP element or other high-copy repeat from passing just because 1-2
    chance bases extend it into a "unique" longer match, while still allowing a
    real att site that has one near-copy elsewhere (seen in the test genomes).
The longest kept candidate is reported (ties: closest to the region ends).
Its outer coordinates give a refined element boundary: dr_left_start..dr_right_end.

Output (in --outdir): dr.tsv, one row per region; versions.txt
"""
import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import reverse_complement

from common import parse_coordinates, read_fasta_list, read_tsv, write_tsv, write_versions


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


TRIMS = [(i, j) for i in range(3) for j in range(3)]   # bp trimmed from (left, right)


def find_all(pattern, text):
    """Start positions of pattern in text, overlapping allowed (fast C-level search)."""
    hits, pos = [], text.find(pattern)
    while pos != -1:
        hits.append(pos)
        pos = text.find(pattern, pos + 1)
    return hits


def trim_counts(pattern, contigs):
    """Copies of every 0-2 bp per-end trim of pattern, on the forward strand.

    Searches only for the core (pattern minus 2 bp per end), then checks the
    bases around each core hit, so one genome scan covers all 9 trims.
    """
    n = len(pattern)
    core = pattern[2:n - 2]
    counts = dict.fromkeys(TRIMS, 0)
    for c in contigs:
        for h in find_all(core, c):
            start = h - 2                       # where the untrimmed pattern would start
            for i, j in TRIMS:
                a, b = start + i, start + n - j
                if a >= 0 and b <= len(c) and c[a:b] == pattern[i:n - j]:
                    counts[(i, j)] += 1
    return counts


def dr_copies(seq, contigs):
    """(copies of seq, max copies of seq trimmed by up to 2 bp per end), both strands."""
    n = len(seq)
    fwd = trim_counts(seq, contigs)
    rev = trim_counts(reverse_complement(seq), contigs)
    total = {}
    for i, j in TRIMS:
        v = seq[i:n - j]
        # trimming i bp left / j bp right of seq = trimming j left / i right of its rc;
        # a reverse-palindromic variant is found on both strands at once: count it once
        total[(i, j)] = fwd[(i, j)] + (0 if reverse_complement(v) == v else rev[(j, i)])
    exact = total[(0, 0)]
    trimmed = max(v for k, v in total.items() if k != (0, 0))
    return exact, trimmed


def count_copies(seq, contigs):
    """Occurrences (overlapping allowed) on both strands of all contigs."""
    rc = reverse_complement(seq)
    strands = [seq] if rc == seq else [seq, rc]
    return sum(len(find_all(x, c)) for c in contigs for x in strands)


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
            seen[s] = dr_copies(s, all_contigs)
        exact, trimmed = seen[s]
        if exact > args.max_copies or trimmed > args.max_copies + 1:
            continue
        # distance of each copy to its region end (for tie-breaking)
        dist = abs(lraw[i] - start) + abs(rraw[j + mlen - 1] - end)
        kept.append((mlen, -dist, i, j, s, exact))
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
    ap.add_argument("--outside", type=int, default=5000,
                    help="bp searched outside each region end (default: 5000)")
    ap.add_argument("--inside", type=int, default=5000,
                    help="bp searched inside each region end (default: 5000)")
    ap.add_argument("--min-len", type=int, default=15,
                    help="minimum DR length (default: 15, as in Ares-Arroyo et al. 2026; "
                         "see README for background rates)")
    ap.add_argument("--max-copies", type=int, default=2,
                    help="max genome-wide copies of a DR (default: 2)")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    regions = read_tsv(args.regions)
    genomes, circular = read_fasta_list(args.fasta_list)

    rows = []
    for genome in sorted({r["genome"] for r in regions}):
        contigs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(genomes[genome], "fasta")}
        all_seqs = list(contigs.values())
        for r in (r for r in regions if r["genome"] == genome):
            segs = parse_coordinates(r["coordinates"])
            start, end = segs[0][0], segs[-1][1]
            if end < start:  # region wraps the origin: express end past the contig length
                end += len(contigs[r["contig"]])
            res = find_dr(contigs[r["contig"]], all_seqs, start, end, args,
                          r["contig"] in circular)
            rows.append({"region_id": r["region_id"], **res})

    write_tsv(out / "dr.tsv", rows, FIELDS)
    write_versions(out / "versions.txt", outside_bp=args.outside, inside_bp=args.inside,
                   min_len=args.min_len, max_copies=args.max_copies,
                   seed="not applicable (deterministic)")
    print(f"DRs found in {sum(r['dr'] for r in rows)} of {len(rows)} regions -> {out / 'dr.tsv'}")


if __name__ == "__main__":
    main()
