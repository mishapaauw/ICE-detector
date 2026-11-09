# Project ICE-detector

This tool delineates ICEs in bacterial genome assemblies. See SPEC.md

The algorithm should be the following:

- Run PPanggolin to make a pangenome and infer 'spots', regions of insertion in particular genomes.  
- For each spot: detect conjugation systems (macsyfinder) and relaxases, oriTs, and GC content.

## Milestones

1. Succesfully running ppanggolin on test genome set, generating spot genomic fasta files.
2. Implementing script that works per spot to detect conjugation systems (2a), relaxases (2b), oriTs (2c), and GC content (2d).
3. Implementing a final output table with high-confidence ICE spots including the amount of evidence you have.

## Layout

- data/test/ read-only inputs, never modify. I will add 5 test genomes here
- results/ generated result tables (gitignored)
- scripts/ scripts to run tools (.sh, or a python script orchestrating everything)
- tests/ i have no experience with this, you'll have to help me here. but it's not the most important part for me.
- environment.yml
- Feel free to extend the layout within reason; adhere to bioinformatics project layouts. 

## Rules

- Don't add dependencies or install anything outside conda without asking.
- Long jobs: run in the background, log to a file, check sparingly.
- Never print long files (e.g. full genome fasta). Use head, wc -l, cut, or write a summary script.
- Set seeds and record tool versions in every output folder.
- Don't change expected test outputs to make tests pass; ask me instead.
- Input genome FASTA files (chromosome + plasmids) are only read, never changed: other tools in the larger pipeline (e.g. geNomad) need them as they are. ICE-detector is standalone (no geNomad output needed): it analyses only the largest contig per genome, via temporary copies that are deleted after the run. See README "Input genomes".

## Environments

- Create a dedicated conda environment named `ICE-detector`; never touch `base`.
- Use conda with channels bioconda and conda-forge
- Tools: ppanggolin, HMMER, and macsyfinder (use pip for this one, if required)
- If tools conflict, use separate environments and tell me why.
- Verify each tool installation
- Export the result to `environment.yml` (pinned versions) and commit it.
- Don't install anything outside conda without asking me.

## Domain specifics

- Use PPanGGOlin v2 to make a pangenome. Use the whole genomes as input (`ppanggolin all --fasta genomes.fasta.list`). You make the fasta list based on the genomes in `data/test/`
