# ICE-detector

<img src="docs/logo.png" alt="ICE-detector logo" width="200">

A tool to detect ICEs in bacterial genomes using pangenomics and function annotation.
Co-developed by Claude Code and myself.

Delineates **integrative conjugative elements (ICEs)** in complete bacterial
genome assemblies of one species. The same evidence is also used to flag
**integrative mobilizable elements (IMEs)**, decayed ICEs and elements that carry
only an oriT. 

Approach (after Cury et al. 2017 and Ares-Arroyo et al. 2026; see [References](#references)):

1. **PPanGGOLiN** builds a pangenome from all genomes and predicts regions of
   genomic plasticity (RGPs), grouped into *spots* (insertion sites shared across genomes).
2. Neighbouring RGPs on the same chromosome that are at most `--merge-gap` (3)
   genes apart are merged into one region: an element is often split into
   several RGPs by a gene that is persistent in the pangenome. RGPs are not
   merged across a tRNA/tmRNA gene, which usually marks the integration site of a
   separate element. A merged region keeps the ID of its largest RGP; `rgp_ids`
   lists all parts.
   Region edges are then **extended over mobility genes**: if one of the 3 genes
   (`--extend-gap`) just outside an edge carries a relaxase, a conjugation-system
   gene (MacSyFinder) or a tyrosine/serine integrase, the edge moves outward to
   include it, and the search repeats from the new edge. Such genes are often
   called persistent by PPanGGOLiN when the element is common in the genome set,
   and are then left outside the RGP. Extension stops at a tRNA/tmRNA and at genes
   of another region. All later steps (oriT, GC, DRs, tRNA, FASTA, final table)
   use the extended edges; `rgp_start`/`rgp_end` keep the PPanGGOLiN edges.
   Each region with at least `--min-orfs` genes (counted after merging) is scored for:
   - a **conjugation system** (T4SS): MacSyFinder + CONJScan, chromosome models
   - a **relaxase** (MOB families, from the same CONJScan run)
   - an **integrase**: tyrosine (Pfam PF00589) or serine (PF00239, PF07508)
     recombinase, hmmsearch with Pfam gathering thresholds. DDE transposases
     (PF00665) are reported but not counted: IS elements carry them everywhere
     (193 hits in the 5 test genomes vs 58 tyrosine/serine recombinases)
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
   - a **tRNA/tmRNA at a region end**: many ICEs integrate into the 3' end of a
     tRNA or tmRNA gene. tRNAs are predicted with ARAGORN (pyaragorn); a region
     counts when a tRNA lies within 1 kb (`--trna-dist`) of either end, or
     overlaps one of the DR copies
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
   | `none` | no | no | none |

   The `_noOriT` categories exist because the oriT profiles only cover known oriT
   families: many real ICEs and IMEs will have an oriT that is not detected.

Only the chromosome is analysed: plasmids are out of scope (identify them
separately, e.g. with geNomad). See [Input genomes](#input-genomes).

## Installation

```bash
conda config --set channel_priority strict   # recommended
conda env create -f environment.yml          # add --solver libmamba if solving is slow
conda activate ICE-detector
```

`environment.yml` pins only the tools the pipeline calls (PPanGGOLiN, MacSyFinder,
HMMER, pyaragorn, Biopython, numpy, pandas) and solves on macOS and Linux.
`envs/ICE-detector.osx-arm64.lock.yml` is the full export of the tested macOS
(arm64) environment, with every package version, for exact reproduction on that platform.

Download the reference data once (stored in `db/`, not tracked by git; see `db/README.md`):

```bash
msf_data install CONJScan==2.1.0 -t db/macsy_models
python scripts/download_orit_hmms.py -o db/orit
python scripts/download_pfam_hmms.py -o db/pfam
```

## Usage

Input: a folder with one FASTA per **complete** genome (≥ 5 genomes of the same
species). See [Input genomes](#input-genomes).

```bash
python scripts/run_all.py data/test -o results/test --cpu 8
```

| Option | Default | Meaning |
|---|---|---|
| `--min-orfs` | 10 | keep only regions with at least this many genes (after merging) |
| `--merge-gap` | 3 | merge neighbouring RGPs at most this many genes apart; 0 = no merging |
| `--extend-gap` | 3 | extend region edges over a relaxase / conjugation / integrase gene at most this many genes outside the edge; 0 = off |
| `--min-evidence` | 2 | minimum evidence count (of 7) for `ice_candidates.tsv` |
| `--orit-min-cov` | 0.7 | an oriT hit must cover this fraction of its HMM profile |
| `--orit-max-evalue` | 0.01 | maximum E-value of an oriT hit |
| `--gc-flank` | 35000 | bp up- and downstream used for flank GC |
| `--gc-z` | 2 | \|GC z-score\| at or above which GC counts as atypical |
| `--dr-outside`, `--dr-inside` | 5000, 5000 | bp searched for DRs outside / inside each region end |
| `--dr-min-len` | 15 | minimum DR length (as in Ares-Arroyo et al. 2026) |
| `--dr-max-copies` | 2 | max genome-wide copies of a DR (versions trimmed by ≤ 2 bp per end may have one more) |
| `--trna-dist` | 1000 | bp between a tRNA/tmRNA and a region end |
| `--linear` | off | treat contigs as linear (default: circular) |
| `--skip-ppanggolin`, `--skip-macsyfinder`, `--skip-nhmmscan`, `--skip-hmmsearch` | off | reuse earlier tool output (e.g. to change only the filters) |

## Input genomes


1. **Input FASTA files are only read**, never modified, renamed or moved.
2. **Only the largest contig of each genome is analysed** (the chromosome, in a
   complete assembly). All other contigs (plasmids, phage contigs, small
   leftovers) are dropped, whatever their FASTA header says. If two contigs are
   equally large, the first one in the file is kept.
3. The largest contig is written to a **temporary copy**,
   `<outdir>/tmp_chromosomes/<genome>.fasta`, which all steps use. The copies are
   deleted when the run ends, also if it fails. Re-runs (e.g. with
   `--skip-ppanggolin`) make them again, identically.
4. `<outdir>/input_contigs.tsv` records for every genome the input file, the kept
   contig and its length, the assembly length and the dropped contigs
   (`name:length`).
5. **Complete genomes only.** If the kept contig is less than 50% of the assembly
   (likely a draft genome), a warning is written to `logs/make_fasta_list.log`;
   the run continues, but most of that genome is then ignored.
6. Genome name = file name without extension; the chromosome names must be unique
   across genomes (plasmid names may repeat). Contigs are treated as circular
   unless `--linear` is given.

`all_regions.tsv` and `genome_gc` describe
chromosomes only, and DR uniqueness is counted on the chromosome only.

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
| `spot_id`, `spot_n_genomes` | PPanGGOLiN spot, and how many genomes have a region in it (spot labels change between runs; compare runs on `region_id`). Merged regions list all their spots (`spot_3;spot_7`) and give the largest genome count |
| `rgp_ids`, `n_merged` | PPanGGOLiN RGPs merged into this region (`;`-joined, genomic order) and how many |
| `rgp_start`, `rgp_end` | region edges as PPanGGOLiN made them (before merging/extension changed `start`/`end`) |
| `extended_left`, `extended_right`, `extended_by` | genes added at each edge by extension, and the mobility genes that triggered it (e.g. `MOBQ`, `Tyr`, `T4SS`) |
| `category` | `ICE`, `dICE`, `IME`, `iOriT`, `ICE_noOriT`, `dICE_noOriT`, `IME_noOriT`, `other`, `none` (see table above) |
| `mobility_genes` | what was found, joined by `+`: `T4SS`, `dCONJ`, `MOB` (relaxase), `oriT`; e.g. `T4SS+MOB` = ICE without a detected oriT |
| `evidence_count` | number of: `t4ss`, `relaxase`, `orit`, `integrase`, `trna`, `gc_atypical`, `dr` |
| `element_start`, `element_end` | outer edges of the DR pair (refined boundaries), if a DR was found |
| `dr_length`, `dr_seq`, `dr_left_offset`, `dr_right_offset` | the DR, and how far each copy lies from the region end (negative = outside the region) |
| `conj_systems` | e.g. `T4SS_typeG:…(21/21)`: system genes inside the region / total |
| `relaxase_types`, `orit_profiles` | MOB families; oriT profiles hit (best E-value) |
| `integrase_types`, `n_integrases`, `integrase_end_dist` | `Tyr` / `Ser` recombinases in the region, how many, and bp from the closest one to a region end |
| `dde_transposase` | DDE transposase (Pfam rve) in the region; information only, not evidence |
| `trna_name`, `trna_end`, `trna_dist`, `trna_dr_overlap` | nearest tRNA/tmRNA (e.g. `tRNA-Gly(gcc)`), which end (`left`/`right`), distance in bp, and whether it overlaps a DR copy |
| `region_gc`, `genome_gc`, `gc_diff`, `gc_zscore` | GC vs the whole genome |
| `flank_up_gc`, `flank_down_gc`, `flank_gc`, `gc_diff_flank` | GC vs the flanking DNA |

Intermediate results are kept per step: `ppanggolin/`, `regions/` (per-region FASTA),
`conjugation/`, `integrase/`, `orit/`, `gc/`, `dr/`, `trna/`, each with its own
`versions.txt`. `trna/trnas.tsv` lists every predicted tRNA/tmRNA. Also in
`<outdir>`: `input_contigs.tsv` (see [Input genomes](#input-genomes)) and
`genomes.fasta.list` (the PPanGGOLiN input; it points to the temporary copies,
which no longer exist after the run).

## Tests

```bash
pytest           # fast unit tests (< 1 s)
pytest -m slow   # full pipeline on data/test vs tests/expected/test_candidates.tsv
```

## Limitations

- Complete genomes only: per genome just the largest contig is analysed, so ICEs
  on plasmids, secondary chromosomes or other contigs are not found, and draft
  assemblies lose most of their sequence.
- Region boundaries come from PPanGGOLiN, refined only when a DR pair is found.
  Neighbouring RGPs are merged by gene distance only, so two separate elements
  that sit within `--merge-gap` genes of each other, with no tRNA between them,
  also end up in one region (check `rgp_ids`).
- Edge extension can pull a conserved system into a small region: a 5-gene RGP
  next to the species-wide type-T T4SS becomes a 15-gene `other` (T4SS) region.
  Check `extended_by` and `rgp_start`/`rgp_end` for such cases. RGPs that span the origin of a circular chromosome are not merged.
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
- Integrase and tRNA evidence are also found in prophages and other islands; like
  GC and DRs they only support a call that has a mobility marker. Integrases of
  the DDE type (some ICE families) are not counted.
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
- Tools: PPanGGOLiN, MacSyFinder with CONJScan models, HMMER, ARAGORN (pyaragorn);
  oriT HMM profiles from figshare doi:10.6084/m9.figshare.31047685; Pfam profiles
  for integrases (see `db/README.md`).

## Layout

```
config/       PPanGGOLiN parameters (seed)
data/test/    test genomes (read-only)
db/           downloaded models (gitignored except README.md)
envs/         full export of the tested conda env
scripts/      pipeline steps + run_all.py (shared helpers in common.py)
tests/        pytest tests and expected outputs
results/      outputs (gitignored)
```
