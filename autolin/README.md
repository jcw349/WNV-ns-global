# Generate WNV lineage names with Autolin

## Background Nomenclature
### Clades and Lineages
West Nile virus have some established nomenclature in previous publications. Main lineages discovered thus far are denoted with roman numerals {i.e. I, II, III, ...} and clades have been identified. Most publications have referred to this and lineage/clade designations {i.e. Ia, Ib, II, ...}

### Strains
There are multiple strains designated throughout various publications. Most commonly found in Lineage-clade Ia, NY99, WN02, SW03 have been used since early 2000. 

### Designations 
To maintain lineage designations with nomenclature that have been used to describe and characterize WNV, we will use lineage-clade.strain to guide lineage names.

## Procedures
### A) Create a simple MATree
#### Option 1: Convert Auspice json to MAT protobuf
_This option is still being explored as [matUtils](https://usher-wiki.readthedocs.io/en/latest/matUtils.html) is still being developed to convert non-COVID auspice trees to MAT protobuf trees. Issue is being tracked: [github.com/yatisht/usher, issues: 378](https://github.com/yatisht/usher/issues/378)_

#### Option 2: Create Simple MAT protobuf

To create MATree, [UShER](https://usher-wiki.readthedocs.io/en/latest/Installation.html) is required. 

At the minimum, a VCF and newick (nwk) tree are required to build a MATree profobuf file.
1. VCF: The auspice snakemake pipeline for the StaPH-B WNV-global build includes a step generating a "results/merged.vcf" file using faToVcf on aligned.fasta and NC_009942.1 WNV reference.
2. nwk tree: The auspice snakemake pipeline for the StaPH-B WNV-global build creates a "results/tree.nwk" file that has been refined using `augur refine`

`usher -t ./results/tree.nwk -v ./results/merged.vcf -o usher_wnv.pb -d ./results/`

### B) Add clade annotations
To add annotations to MATree, [matUtils annotate](https://usher-wiki.readthedocs.io/en/latest/matUtils.html#annotate) is required.

If a list of lineage designations are already known, use the "config/clades.tsv" file to guide the autolin lineage designation. The file format is tab-delimited with `clades` {i.e. Ia.4.WN02} in the first column and the `strain` in the second column {i.e. KJ501251} (this is the first column of the augur metadata).
Otherwise skip this step to start the lineage designation from A.0, A.1, A.2, etc.

`matUtils annotate -i ./results/usher_wnv.pb -o ./results/usher_wnv_clades.pb -c ./config/clades.tsv` 

### C) Run autolin to designate lineages
Once the clade designations have been annotated onto the MATree, [autolin](https://github.com/jmcbroome/autolin) can now use the annotated clades in the MATree to generate new labels for sublineages.

Follow installation instructions to create conda environment. 

`conda activate autolin`

`wget -P ./scripts https://github.com/jmcbroome/autolin/blob/32d9a52a744d762e243cd83449972aab4d2fa556/propose_sublineages.py`

`python3 ./scripts/propose_sublineages.py -i ./results/usher_wnv_clades.pb -o ./results/autolin_wnv.pb -l ./results/autolin_wnv.labels.txt --recursive -m 3 -t 0 -f 0`

A tab-delimited file, "autolin_wnv.labels.txt" will be created, which can be added to the ./results/metadata.tsv under the `lineages` column prior to the final `augur export v2` to create the auspice json file.
