"""Small helpers shared by the pipeline scripts (TSV I/O, fasta list, coordinates)."""
import csv
import sys
from collections import defaultdict
from datetime import datetime, timezone


def read_tsv(path, comment=None):
    """Rows of a TSV with a header, as dicts. Lines starting with `comment` are skipped."""
    with open(path) as fh:
        lines = [l for l in fh if l.strip() and not (comment and l.startswith(comment))]
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    """Write dicts as TSV. Fields default to the keys of the first row."""
    fields = fields or (list(rows[0]) if rows else [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)


def read_fasta_list(path):
    """PPanGGOLiN fasta list -> ({genome: fasta path}, {circular contig names})."""
    genomes, circular = {}, set()
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row:
                genomes[row[0]] = row[1]
                circular.update(row[2:])
    return genomes, circular


def parse_coordinates(coords):
    """'100..200' or '5000..5100,1..50' (wraps the origin) -> [(100, 200), ...]."""
    return [tuple(map(int, part.split(".."))) for part in coords.split(",")]


def region_segments(regions):
    """contig -> [(start, end, region_id)], one entry per segment of each region."""
    segs = defaultdict(list)
    for r in regions:
        for s, e in parse_coordinates(r["coordinates"]):
            segs[r["contig"]].append((s, e, r["region_id"]))
    return segs


def find_region(segs, contig, start, end):
    """ID of the first region overlapping start..end on contig, or ''."""
    for s, e, rid in segs.get(contig, []):
        if start <= e and end >= s:
            return rid
    return ""


def write_versions(path, **settings):
    """versions.txt: date, python version, then `key: value` per setting."""
    with open(path, "w") as fh:
        fh.write(f"date: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n"
                 f"python: {sys.version.split()[0]}\n")
        for k, v in settings.items():
            fh.write(f"{k}: {v}\n")
