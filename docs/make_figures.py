#!/usr/bin/env python
"""Draw docs/workflow.png and docs/logo.png (matplotlib only).

    python docs/make_figures.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Wedge

OUT = Path(__file__).resolve().parent

INK = "#1f2933"
MUTED = "#52606d"
ICE = "#2b8cbe"       # element / ICE colour
ICE_LIGHT = "#d6ecf7"
STEP = "#f0f4f8"
EVID = "#fff4e0"
EVID_EDGE = "#d9822b"
OUTC = "#e3f4e8"
OUT_EDGE = "#3f9142"


def box(ax, x, y, w, h, title, body="", face=STEP, edge=MUTED, title_size=11, body_size=8.6):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=face, ec=edge, lw=1.4))
    if body:
        ax.text(x, y + h / 2 - 0.2, title, ha="center", va="top", fontsize=title_size,
                weight="bold", color=INK)
        ax.text(x, y + h / 2 - 0.55, body, ha="center", va="top", fontsize=body_size,
                color=MUTED, linespacing=1.35)
    else:
        ax.text(x, y, title, ha="center", va="center", fontsize=title_size, weight="bold",
                color=INK)


def arrow(ax, x1, y1, x2, y2, color=MUTED):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.4, color=color, shrinkA=0, shrinkB=0))


def workflow():
    fig, ax = plt.subplots(figsize=(13, 11.5))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 11.5)
    ax.axis("off")

    ax.text(6.5, 11.2, "ICE-detector workflow", ha="center", va="center",
            fontsize=18, weight="bold", color=INK)
    ax.text(6.5, 10.8, "delineating ICEs and IMEs in complete genomes of one species",
            ha="center", va="center", fontsize=10.5, color=MUTED)

    cx = 6.5
    box(ax, cx, 9.95, 7.2, 1.0, "Input: ≥ 5 complete genomes (FASTA)",
        "chromosome + plasmids, one file per genome; files are only read, never changed")
    arrow(ax, cx, 9.45, cx, 9.15)
    box(ax, cx, 8.6, 7.2, 1.1, "1 · Keep the chromosome",
        "largest contig per genome → temporary copies (deleted after the run)\n"
        "kept / dropped contigs recorded in input_contigs.tsv")
    arrow(ax, cx, 8.05, cx, 7.75)
    box(ax, cx, 7.2, 7.2, 1.1, "2 · Pangenome  (PPanGGOLiN)",
        "regions of genomic plasticity (RGPs), grouped into spots (shared insertion sites)\n"
        "seeded (42), versions recorded in every output folder")
    arrow(ax, cx, 6.65, cx, 6.35)
    box(ax, cx, 5.95, 7.2, 0.8, "3 · Regions with ≥ 10 genes  →  sequence + proteins", "")

    # evidence boxes
    ev = [
        ("Conjugation", "system: MacSyFinder\n+ CONJScan\nT4SS (complete)\ndCONJ (decayed)"),
        ("Relaxase", "MOB families\n(same CONJScan run)"),
        ("oriT", "nhmmscan,\n129 oriT HMMs\nE ≤ 0.01,\n≥ 70% profile"),
        ("GC content", "z-score vs genome\n|z| ≥ 2 = atypical\n+ 35 kb flanks"),
        ("Direct repeats", "attL / attR\n≥ 15 bp, ±5 kb\naround the ends,\nunique in genome"),
    ]
    w, gap = 2.3, 0.18
    xs = [cx + (i - 2) * (w + gap) for i in range(5)]
    ax.add_patch(FancyBboxPatch((xs[0] - w / 2 - 0.12, 3.12), 5 * w + 4 * gap + 0.24, 2.18,
                                boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc="none", ec=EVID_EDGE, lw=1, ls=(0, (4, 3))))
    ax.text(xs[0] - w / 2, 5.18, "4 · evidence per region", fontsize=9, color=EVID_EDGE,
            weight="bold", va="center")
    for x, (t, b) in zip(xs, ev):
        arrow(ax, cx, 5.55, x, 5.05, color="#9aa5b1")
        box(ax, x, 4.1, w, 1.7, t, b, face=EVID, edge=EVID_EDGE, title_size=10.5,
            body_size=8.4)
    for x in xs:
        arrow(ax, x, 3.25, cx, 2.55, color="#9aa5b1")

    box(ax, cx, 2.0, 7.2, 1.1, "5 · Combine",
        "evidence count (0-5)  ·  category from oriT, relaxase and conjugation system\n"
        "candidate = evidence ≥ 2 and at least one mobility gene (T4SS, MOB or oriT)")
    arrow(ax, cx, 1.45, cx, 1.15)
    box(ax, cx, 0.68, 7.2, 0.9, "Output",
        "ice_candidates.tsv  ·  all_regions.tsv  ·  ice_fasta/  ·  versions.txt",
        face=OUTC, edge=OUT_EDGE)

    # category legend (right) — kept clear of the evidence boxes
    lx, ly = 10.45, 10.25
    ax.text(lx, ly, "Categories", fontsize=10, weight="bold", color=INK)
    rows = [("", "oriT", "MOB", "T4SS"),
            ("ICE", "+", "+", "complete"),
            ("dICE", "+", "+", "decayed"),
            ("IME", "+", "+", "–"),
            ("iOriT", "+", "–", "–"),
            ("…_noOriT", "–", "+", "as above"),
            ("other", "any other", "", "")]
    cols = [lx, lx + 0.95, lx + 1.45, lx + 1.95]
    for i, r in enumerate(rows):
        y = ly - 0.32 - i * 0.27
        for j, (c, v) in enumerate(zip(cols, r)):
            ax.text(c, y, v, fontsize=8, color=INK if i else MUTED,
                    weight="bold" if (j == 0 or i == 0) else "normal", va="center")

    # mini element sketch (left)
    sx, sy = 0.35, 9.45
    ax.text(sx, sy + 0.8, "What is searched for", fontsize=10, weight="bold", color=INK)
    ax.plot([sx, sx + 2.4], [sy, sy], color="#9aa5b1", lw=4, solid_capstyle="butt")
    ax.plot([sx + 0.6, sx + 1.8], [sy, sy], color=ICE, lw=8, solid_capstyle="butt")
    for x in (sx + 0.6, sx + 1.8):
        ax.add_patch(Polygon([[x - 0.09, sy + 0.12], [x + 0.09, sy + 0.12], [x, sy + 0.27]],
                             color=EVID_EDGE))
    ax.text(sx + 1.2, sy - 0.42, "ICE / IME", ha="center", va="top", fontsize=8.5, color=ICE,
            weight="bold")
    ax.text(sx + 0.6, sy + 0.33, "attL", ha="center", fontsize=7.5, color=EVID_EDGE)
    ax.text(sx + 1.8, sy + 0.33, "attR", ha="center", fontsize=7.5, color=EVID_EDGE)
    ax.text(sx, sy - 0.12, "chromosome", va="top", fontsize=7.5, color=MUTED)

    fig.savefig(OUT / "workflow.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def logo():
    fig, ax = plt.subplots(figsize=(4, 4.6))
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.75, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    # circular chromosome with an inserted element (arc) between two direct repeats
    ax.add_patch(Circle((0, 0), 0.85, fill=False, ec="#9aa5b1", lw=11))
    start, end = 20, 110
    ax.add_patch(Arc((0, 0), 1.7, 1.7, theta1=start, theta2=end, color=ICE, lw=15))
    for a in (start, end):
        t = np.radians(a)
        ax.add_patch(Wedge((0.85 * np.cos(t), 0.85 * np.sin(t)), 0.13, 0, 360,
                           color=EVID_EDGE))
    # magnifying glass over the element
    gx, gy = 0.18, 0.2
    ax.add_patch(Circle((gx, gy), 0.42, fc=ICE_LIGHT, ec=INK, lw=6, alpha=0.95))
    ax.plot([gx - 0.3, gx - 0.75], [gy - 0.3, gy - 0.75], color=INK, lw=12,
            solid_capstyle="round")
    # snowflake-like sparkle inside the lens (the "ice")
    for k in range(3):
        a = np.radians(90 + 60 * k)
        dx, dy = 0.24 * np.cos(a), 0.24 * np.sin(a)
        ax.plot([gx - dx, gx + dx], [gy - dy, gy + dy], color=ICE, lw=4, solid_capstyle="round")
    ax.add_patch(Circle((gx, gy), 0.05, color=ICE))

    # centre the two-coloured name: measure both parts, then place them
    kw = dict(va="center", fontsize=26, weight="bold")
    a = ax.text(0, -1.38, "ICE", color=ICE, **kw)
    b = ax.text(0, -1.38, "-detector", color=INK, **kw)
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    width = lambda t: np.diff(inv.transform(t.get_window_extent(r))[:, 0])[0]
    wa, wb = width(a), width(b)
    a.set_x(-(wa + wb) / 2)
    b.set_x(-(wa + wb) / 2 + wa)
    fig.savefig(OUT / "logo.png", dpi=150, bbox_inches="tight", transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    workflow()
    logo()
    print(f"-> {OUT / 'workflow.png'}, {OUT / 'logo.png'}")
