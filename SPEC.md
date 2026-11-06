## Goal of this tool

The goal of this tool is to delineate integrative conjugative elements in bacterial genome assemblies. It will use a 2-fold strategy: delineating spots of genome flexibility in a bacterial pangenome using PPanggolin, and functionally annotation them for the presence of a conjugation system and a oriT.

## Inspiration

From Cury, Touchon, Rocha (2017):

> "ICEs were delimited in two steps. First, we identified the spots encoding conjugative systems (see definition above). The core genes flanking these spots provide upper bounds for the limits of the ICE. Second, we analyzed the interval with the ICE and identified the limits of the element by overlaying the information on the presence of genes of the conjugation system, on G+C content, and on the frequency of accessory genes in the spot pan-genome. The genes of the ICE are expected to be present in the spot pan-genome at similar frequencies (some differences may be caused by mutations, deletions, transposable elements, and annotation errors) and this information is usually sufficient to delimit the ICE."

## Input:

A folder with bacterial genome assemblies of the same species (minimum amount: 5, maximum: ?)

## Output

- A table with genomic spots: genome ID, contig ID, start and end coordinates, number of ORFs, added with T4SS presence, oriT presence, relaxase presence.
- A folder with the genomic sequences of the spots (taken from the input genomes). .fasta format 

## Options

- There should be an option to only output spots made at least <n> ORFs. This is useful to not save many small 1 or 2 gene insertions.

## Tools

- Spot detection: PPanGGolin
- Conjugation systems detection: Macsyfinder with CONJscan models
- OriT detection: HMMer search using oriT profiles from figshare (https://figshare.com/articles/dataset/HMM_profiles_for_the_screening_of_conjugative_origins_of_transfer/31047685), first HMMpress.

## References

> Cury, J., Touchon, M., Rocha, E.. 2017, Nucleic Acid Research. doi: https://doi.org/10.1093/nar/gkx607
