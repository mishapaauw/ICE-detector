# Reference databases (not tracked in git)

Everything in `db/` except this file is downloaded and gitignored.

## CONJScan (MacSyFinder models) — `db/macsy_models/CONJScan`
- Version: 2.1.0
- Source: MacSyFinder model registry (github.com/macsy-models/CONJScan)
- Install (inside the `ICE-detector` env):
  `msf_data install CONJScan==2.1.0 -t db/macsy_models`
- Used models: `CONJScan/Chromosome all` (T4SS_*, dCONJ_*, MOB)

## oriT HMM profiles — `db/orit/` 
- Source: https://figshare.com/articles/dataset/HMM_profiles_for_the_screening_of_conjugative_origins_of_transfer/31047685
- Version: figshare v1 (doi 10.6084/m9.figshare.31047685.v1), 129 DNA profiles
- Download + MD5 check + concatenate: `python scripts/download_orit_hmms.py -o db/orit`
- `detect_orit.sh` runs `hmmpress` on `db/orit/oriT.hmm` if needed

## Pfam integrase/recombinase profiles — `db/pfam/`
- Source: InterPro API (https://www.ebi.ac.uk/interpro/api/), Pfam release recorded
  in `db/pfam/manifest.tsv` (39.0 when first downloaded, 2026-11-11)
- Profiles: PF00589 Phage_integrase (Tyr), PF00239 Resolvase (Ser),
  PF07508 Recombinase (Ser), PF00665 rve (DDE, reported only)
- Download + concatenate + hmmpress: `python scripts/download_pfam_hmms.py -o db/pfam`
- Searched with `hmmsearch --cut_ga` (Pfam gathering thresholds)
