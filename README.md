# ICE-detector

A tool to detect ICEs in bacterial genomes using pangenomics and function annotation.
Co-developed by Claude Code and myself.

Delineates **integrative conjugative elements (ICEs)** in complete or draft bacterial
genome assemblies of one species. The same evidence is also used to flag
**integrative mobilizable elements (IMEs)**, decayed ICEs and elements that carry
only an oriT. Finding these is an explicit goal too, although the tool keeps its
name for now.

Approach (after Cury et al. 2017 and Ares-Arroyo et al. 2026; see [References](#references)):

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
     region end is searched 5 kb outward and 5 kb inward for an exact,
     same-strand repeat of ≥ 15 bp that occurs only at these two places in the
     genome. When found, the DR pair gives a refined element boundary
     (`element_start`/`element_end`)
3. Evidence is combined into a table with a category (following Ares-Arroyo et
   al. 2026) and an evidence count:

   | Category | oriT | Relaxase | Conjugation system |
   |---|---|---|---|
   | `ICE` | yes | yes | complete (CONJScan T4SS model) |
   | `dICE` (decayed ICE) | yes | yes | decayed (CONJScan dCONJ model) |
   | `IME` | yes | yes | none |
   | `iOriT` | yes | no | none |
   | `ICE_noOriT` | not found | yes | complete |
   | `dICE_noOriT` | not found | yes | decayed |
   | `IME_noOriT` | not found | yes | none |
   | `other` | any other combination (e.g. conjugation system without relaxase); see `mobility_genes` | | |

   The `_noOriT` categories exist because the oriT profiles only cover known oriT
   families: many real ICEs and IMEs will have an oriT that is not detected.
   | `none` | no | no | none |

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
| `--dr-outside`, `--dr-inside` | 5000, 5000 | bp searched for DRs outside / inside each region end |
| `--dr-min-len` | 15 | minimum DR length (as in Ares-Arroyo et al. 2026) |
| `--dr-max-copies` | 2 | max genome-wide copies of a DR (versions trimmed by ≤ 2 bp per end may have one more) |
| `--linear` | off | treat contigs as linear (default: circular) |
| `--skip-ppanggolin`, `--skip-macsyfinder`, `--skip-nhmmscan` | off | reuse earlier tool output (e.g. to change only the filters) |

## Output (`<outdir>/final/`)

- `ice_candidates.tsv`: regions with evidence count ≥ `--min-evidence` **and** at least one
  mobility marker (T4SS, relaxase or oriT): ICE, dICE, IME and iOriT candidates, plus
  `other` combinations
- `all_regions.tsv`: every region with all evidence columns
- `ice_fasta/`: genomic sequence of each candidate region
- `versions.txt`: tool versions, seeds and settings of every step

Key columns:

| Column | Meaning |
|---|---|
| `genome`, `contig`, `start`, `end`, `length`, `n_orfs` | location of the region |
| `spot_id`, `spot_n_genomes` | PPanGGOLiN spot, and how many genomes have a region in it (spot labels change between runs; compare runs on `region_id`) |
| `category` | `ICE`, `dICE`, `IME`, `iOriT`, `ICE_noOriT`, `dICE_noOriT`, `IME_noOriT`, `other`, `none` (see table above) |
| `mobility_genes` | what was found, joined by `+`: `T4SS`, `dCONJ`, `MOB` (relaxase), `oriT`; e.g. `T4SS+MOB` = ICE without a detected oriT |
| `evidence_count` | number of: `t4ss`, `relaxase`, `orit`, `gc_atypical`, `dr` |
| `element_start`, `element_end` | outer edges of the DR pair (refined boundaries), if a DR was found |
| `dr_length`, `dr_seq`, `dr_left_offset`, `dr_right_offset` | the DR, and how far each copy lies from the region end (negative = outside the region) |
| `conj_systems` | e.g. `T4SS_typeG:…(21/21)`: system genes inside the region / total |
| `relaxase_types`, `orit_profiles` | MOB families; oriT profiles hit (best E-value) |
| `region_gc`, `genome_gc`, `gc_diff`, `gc_zscore` | GC vs the whole genome |
| `flank_up_gc`, `flank_down_gc`, `flank_gc`, `gc_diff_flank` | GC vs the flanking DNA |

Intermediate results are kept per step: `ppanggolin/`, `regions/` (per-region FASTA),
`conjugation/`, `orit/`, `gc/`, `dr/`, each with its own `versions.txt`.

## Tests

```bash
pytest           # fast unit tests (< 1 s)
pytest -m slow   # full pipeline on data/test vs tests/expected/test_candidates.tsv
```

## Limitations

- Region boundaries come from PPanGGOLiN, refined only when a DR pair is found.
  Refinement using gene frequencies within the spot (as in Cury et al. 2017) is
  not implemented.
- DRs: only exact repeats are found. A diverged or deleted attR copy (seen in one
  test ICE) gives no DR. Short DRs also occur by chance. On the test genomes,
  fake regions placed at random in the conserved genome got a "DR" in:

  | min length | 3 kb windows | 5 kb windows (default) |
  |---|---|---|
  | 15 bp (default) | 6.7% | 12.0% |
  | 16 bp | 5.7% | 11.0% |
  | 18 bp | 4.0% | 6.7% |
  | 20 bp | 2.0% | 3.7% |

  vs 40% of the real regions. Treat a DR as supporting evidence, and prefer DRs
  whose copies lie close to the region ends (`dr_left_offset`, `dr_right_offset`).
- DR + atypical GC also occurs in elements without mobility genes (e.g. prophages
  or other genomic islands). These are kept out of `ice_candidates.tsv` (a
  mobility marker is required) but remain in `all_regions.tsv`.
- MacSyFinder can attach a relaxase far from its T4SS (relaxases are "loners" in
  CONJScan); the `(n/total)` counts in `conj_systems` show such partial overlaps.

## References

- Cury J, Touchon M, Rocha EPC (2017). Integrative and conjugative elements and
  their hosts: composition, distribution and organization. *Nucleic Acids Research*
  45(15):8943–8956. doi:10.1093/nar/gkx607
- Ares-Arroyo M, Junker R, Nucci A, Touchon M, Rocha EPC (2026). Genomic Islands as
  minimal hitchers of conjugative elements. *bioRxiv* (v2).
  doi:10.64898/2026.01.13.699239. Source of the ICE / dICE / IME / iOriT categories
  and of the 15 bp minimum DR length.
- Tools: PPanGGOLiN, MacSyFinder with CONJScan models, HMMER; oriT HMM profiles from
  figshare doi:10.6084/m9.figshare.31047685 (see `db/README.md`).

## Layout

```
config/       PPanGGOLiN parameters (seed)
data/test/    test genomes (read-only)
db/           downloaded models (gitignored except README.md)
scripts/      pipeline steps + run_all.py
tests/        pytest tests and expected outputs
results/      outputs (gitignored)
```
