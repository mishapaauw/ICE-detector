# Reference databases (not tracked in git)

Everything in `db/` except this file is downloaded and gitignored.

## CONJScan (MacSyFinder models) — `db/macsy_models/CONJScan`
- Version: 2.1.0
- Source: MacSyFinder model registry (github.com/macsy-models/CONJScan)
- Install (inside the `ICE-detector` env):
  `msf_data install CONJScan==2.1.0 -t db/macsy_models`
- Used models: `CONJScan/Chromosome all` (T4SS_*, dCONJ_*, MOB)

## oriT HMM profiles — `db/orit/` (milestone 4)
- Source: https://figshare.com/articles/dataset/HMM_profiles_for_the_screening_of_conjugative_origins_of_transfer/31047685
