#!/usr/bin/env python
"""Download the oriT HMM profiles from figshare, check MD5s, write one HMM db.

Source: https://figshare.com/articles/dataset/31047685 (doi 10.6084/m9.figshare.31047685)
Output (in --outdir):
  raw/<name>          downloaded files, as published
  oriT.hmm            all profiles concatenated (run `hmmpress` on it)
  manifest.tsv        file name, size, md5, figshare url
"""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

ARTICLE = 31047685


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "ICE-detector"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="db/orit")
    ap.add_argument("--version", type=int, default=1, help="figshare article version")
    args = ap.parse_args()

    out = Path(args.outdir)
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    meta = json.loads(fetch(f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/{args.version}"))

    rows = []
    for f in meta["files"]:
        dest = raw / Path(f["name"]).name  # never trust paths from the remote
        if not (dest.exists() and hashlib.md5(dest.read_bytes()).hexdigest() == f["computed_md5"]):
            data = fetch(f["download_url"])
            if hashlib.md5(data).hexdigest() != f["computed_md5"]:
                raise SystemExit(f"MD5 mismatch for {f['name']}")
            dest.write_bytes(data)
        rows.append((dest.name, f["size"], f["computed_md5"], f["download_url"]))

    hmms = sorted(p for p in raw.iterdir() if p.read_bytes()[:6] == b"HMMER3")
    with open(out / "oriT.hmm", "wb") as fh:
        for p in hmms:
            data = p.read_bytes()
            fh.write(data if data.endswith(b"\n") else data + b"\n")
    with open(out / "manifest.tsv", "w") as fh:
        fh.write(f"# figshare article {ARTICLE} version {args.version} doi {meta.get('doi')}\n")
        fh.write("file\tsize\tmd5\turl\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} files downloaded, {len(hmms)} HMM files -> {out / 'oriT.hmm'}")


if __name__ == "__main__":
    main()
