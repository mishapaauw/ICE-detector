# ICE-detector

A tool to detect ICEs in bacterial genomes using pangenomics and function annotation.
Co-developed by Claude Code and myself.

Delineates **integrative conjugative elements (ICEs)** in complete or draft bacterial
genome assemblies of one species. The same evidence is also used to flag
**integrative mobilizable elements (IMEs)**: elements that carry a relaxase and/or
oriT but no conjugation system of their own. Finding IMEs is an explicit goal too,
although the tool keeps its name for now.

Approach (after Cury, Touchon & Rocha 2017, *NAR*, doi:10.1093/nar/gkx607):

1. **PPanGGOLiN** builds a pangenome from all genomes and predicts regions of
   genomic plasticity (RGPs), grouped into *spots* (insertion sites shared across genomes).
2. Each region with at least `--min-orfs` genes is scored for:
   - a **conjugation system** (T4SS): MacSyFinder + CONJScan, chromosome models
   - a **relaxase** (MOB families, from the same CONJScan run)
   - an **origin of transfer (oriT)**: nhmmscan with oriT DNA profiles, filtered on
     E-value and on how much of the profile the hit covers
   - **GC content**: against the genome (z-score vs same-length windows) and
     against the 35 kb flanking the region on each side
   - **direct repeats (DRs)**: integration duplicates the attachment site, so an
     integrated element is flanked by two identical copies (attL/attR). Each
     region end is searched 3 kb outward and 3 kb inward for an exact,
     same-strand repeat of ≥ 18 bp that occurs only at these two places in the
     genome. When found, the DR pair gives a refined element boundary
     (`element_start`/`element_end`)
3. Evidence is combined into a table with a category (ICE / decayed_ICE /
   mobilizable / none), an evidence count and, for ICEs, a confidence level.

Plasmids are out of scope; identify them separately (e.g. with geNomad).

## Installation

```bash
conda env create -n ICE-detector -f environment.yml
conda activate ICE-detector
```

Download the reference data once (stored in `db/`, not tracked by git; see `db/README.md`):

```bash
msf_data install CONJScan==2.1.0 -t db/macsy_models
python scripts/download_orit_hmms.py -o db/orit
```

## Usage

Input: a folder with one FASTA per genome (≥ 5 genomes of the same species).
Contig names must be unique across all genomes.

```bash
python scripts/run_all.py data/test -o results/test --cpu 8
```

| Option | Default | Meaning |
|---|---|---|
| `--min-orfs` | 10 | keep only regions with at least this many genes |
| `--min-evidence` | 2 | minimum evidence count (of 5) for `ice_candidates.tsv` |
| `--orit-min-cov` | 0.7 | an oriT hit must cover this fraction of its HMM profile |
| `--orit-max-evalue` | 0.01 | maximum E-value of an oriT hit |
| `--gc-flank` | 35000 | bp up- and downstream used for flank GC |
| `--gc-z` | 2 | \|GC z-score\| at or above which GC counts as atypical |
| `--dr-outside`, `--dr-inside` | 3000, 3000 | bp searched for DRs outside / inside each region end |
| `--dr-min-len` | 18 | minimum DR length |
| `--dr-max-copies` | 2 | max genome-wide copies of a DR (median over its 18-mers) |
| `--linear` | off | treat contigs as linear (default: circular) |
| `--skip-ppanggolin`, `--skip-macsyfinder`, `--skip-nhmmscan` | off | reuse earlier tool output (e.g. to change only the filters) |

## Output (`<outdir>/final/`)

- `ice_candidates.tsv`: regions with evidence count ≥ `--min-evidence` **and** at least one
  mobility marker (T4SS, relaxase or oriT): ICEs and IME candidates
- `all_regions.tsv`: every region with all evidence columns
- `ice_fasta/`: genomic sequence of each candidate region
- `versions.txt`: tool versions, seeds and settings of every step

Key columns:

| Column | Meaning |
|---|---|
| `genome`, `contig`, `start`, `end`, `length`, `n_orfs` | location of the region |
| `spot_id`, `spot_n_genomes` | PPanGGOLiN spot, and how many genomes have a region in it (spot labels change between runs; compare runs on `region_id`) |
| `category` | `ICE` (T4SS), `decayed_ICE` (degraded conjugation system), `mobilizable` (relaxase without conjugation system: IME candidate), `none` |
| `confidence` | ICEs only: `high` = T4SS + relaxase + oriT, `medium` = T4SS + one, `low` = T4SS only |
| `evidence_count` | number of: `t4ss`, `relaxase`, `orit`, `gc_atypical`, `dr` |
| `element_start`, `element_end` | outer edges of the DR pair (refined boundaries), if a DR was found |
| `dr_length`, `dr_seq`, `dr_left_offset`, `dr_right_offset` | the DR, and how far each copy lies from the region end (negative = outside the region) |
| `conj_systems` | e.g. `T4SS_typeG:…(21/21)`: system genes inside the region / total |
| `relaxase_types`, `orit_profiles` | MOB families; oriT profiles hit (best E-value) |
| `region_gc`, `genome_gc`, `gc_diff`, `gc_zscore` | GC vs the whole genome |
| `flank_up_gc`, `flank_down_gc`, `flank_gc`, `gc_diff_flank` | GC vs the flanking DNA |

Intermediate results are kept per step: `ppanggolin/`, `regions/` (per-region FASTA),
`conjugation/`, `orit/`, `gc/`, each with its own `versions.txt`.

## Tests

```bash
pytest           # fast unit tests (< 1 s)
pytest -m slow   # full pipeline on data/test vs tests/expected/test_candidates.tsv (~2 min)
```

## Limitations

- Region boundaries come from PPanGGOLiN, refined only when a DR pair is found.
  Refinement using gene frequencies within the spot (as in Cury et al. 2017) is
  not implemented.
- DRs: only exact repeats are found. A diverged or deleted attR copy (seen in one
  test ICE) gives no DR. Short DRs can occur by chance: on the test genomes,
  random same-length windows in the conserved genome gave a "DR" in ~3.5% of
  cases at 18 bp (vs ~29% at 15 bp), which is why 18 bp is the default.
- DR + atypical GC also occurs in elements without mobility genes (e.g. prophages
  or other genomic islands). These are kept out of `ice_candidates.tsv` (a
  mobility marker is required) but remain in `all_regions.tsv`.
- MacSyFinder can attach a relaxase far from its T4SS (relaxases are "loners" in
  CONJScan); the `(n/total)` counts in `conj_systems` show such partial overlaps.

## Layout

```
config/       PPanGGOLiN parameters (seed)
data/test/    test genomes (read-only)
db/           downloaded models (gitignored except README.md)
scripts/      pipeline steps + run_all.py
tests/        pytest tests and expected outputs
results/      outputs (gitignored)
```
