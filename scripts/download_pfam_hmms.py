#!/usr/bin/env python
"""Download the Pfam HMMs used for integrase detection and build one HMM db.

Source: InterPro API (https://www.ebi.ac.uk/interpro/api/), one HMM per Pfam family.
Profiles (role used by summarize_integrase.py):
  PF00589  Phage_integrase  tyrosine recombinase   -> integrase evidence (Tyr)
  PF00239  Resolvase        serine recombinase     -> integrase evidence (Ser)
  PF07508  Recombinase      large serine recomb.   -> integrase evidence (Ser)
  PF00665  rve              DDE transposase        -> reported only (IS elements)
Output (in --outdir):
  raw/<acc>.hmm   downloaded profiles
  integrase.hmm   all profiles concatenated, hmmpress'ed
  manifest.tsv    Pfam release, accession.version, name, length, GA cutoff, md5
"""
import argparse
import gzip
import hashlib
import json
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://www.ebi.ac.uk/interpro/api"
from summarize_integrase import ROLES as PROFILES  # accession -> role


def fetch(url, tries=5):
    """GET with retries; the API answers 429 (too many requests) when called quickly."""
    req = urllib.request.Request(url, headers={"User-Agent": "ICE-detector"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code not in (408, 429, 500, 502, 503, 504) or i == tries - 1:
                raise
            time.sleep(int(e.headers.get("Retry-After") or 0) or 10 * (i + 1))


def header_field(hmm_text, key):
    for line in hmm_text.splitlines():
        if line.startswith(key + " "):
            return line[len(key):].strip()
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="db/pfam")
    args = ap.parse_args()

    out = Path(args.outdir)
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    release = json.loads(fetch(f"{API}/"))["databases"]["pfam"]["version"]

    rows, texts = [], []
    for acc, role in PROFILES.items():
        time.sleep(2)
        data = fetch(f"{API}/entry/pfam/{acc}?annotation=hmm")
        text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
        if not text.startswith("HMMER3") or header_field(text, "ACC").split(".")[0] != acc:
            raise SystemExit(f"{acc}: download is not the expected HMM")
        if not header_field(text, "GA"):
            raise SystemExit(f"{acc}: no GA (gathering) cutoff in the HMM")
        (raw / f"{acc}.hmm").write_text(text)
        texts.append(text if text.endswith("\n") else text + "\n")
        rows.append((header_field(text, "ACC"), header_field(text, "NAME"), role,
                     header_field(text, "LENG"), header_field(text, "GA").rstrip(";"),
                     hashlib.md5(text.encode()).hexdigest()))

    hmm = out / "integrase.hmm"
    hmm.write_text("".join(texts))
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        Path(str(hmm) + ext).unlink(missing_ok=True)
    subprocess.run(["hmmpress", str(hmm)], check=True, stdout=subprocess.DEVNULL)

    with open(out / "manifest.tsv", "w") as fh:
        fh.write(f"# Pfam release {release} via {API}\n")
        fh.write("accession\tname\trole\tlength\tga_cutoff\tmd5\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"Pfam {release}: {len(rows)} profiles -> {hmm}")


if __name__ == "__main__":
    main()
