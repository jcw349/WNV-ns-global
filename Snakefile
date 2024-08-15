rule all:
    input:
        "results/tree.nwk",
        "auspice/tree.json",
        "results/tree.nwk"

#rule usher:
#    input:
#        "results/MATree_annotated.pb",
#        "results/MATree_annot_mutations.tsv",
#        "results/tree.nwk"

#rule autolin:
#    input:
#        "results/MATree_lineages.pb",
#        "results/autolin_lineages.tsv",
#        "results/autolin_labels.tsv"


rule files:
    params:
#        run_id = $(date +"%Y%m%d_%H%M%S"),
        input_fasta = "data/full_dataset.fasta",
        input_metadata = "data/headers.csv",
        reference = "config/reference.gb",
        root = "config/refernece.fasta",
        auspice_config = "config/auspice_config_v2.json",
        clades = "config/clade.tsv",
        known_clades = "config/known_clades.tsv",
        lat_longs = "config/lat_longs.tsv",
        exclude= "config/exclude.txt",
        include= "config/include.txt",
        usher_mutations = "config/nt_clades.tsv",
		lineage_mutations = "config/nt_lineages.tsv"
		tree = "results/tree.nwk"

file = rules.files.params

rule parse:
    message:
        "Parsing {input.sequences}, {input.metadata} and forming FASTA + metadata TSV"
    input:
        sequences = {file.input_fasta},
        metadata = {file.input_metadata}
    output:
        sequences = "results/sequences_raw.fasta",
        metadata = "results/metadata_sans_authors.tsv"
    log:
        "logs/00a-parse-{file.run_id}.log"
    shell:
        """
        python ./scripts/parse_fasta_csv.py \
            {input.sequences} \
            {input.metadata} \
            {output.sequences} \
            {output.metadata} 2> {log}
        """

rule add_authors:
    message:
        "Adding authors to {input.metadata} -> {output.metadata} by collecting info from ENTREZ"
    input:
        metadata = rules.parse.output.metadata
    output:
        metadata = "results/metadata_raw.tsv"
    log:
        "logs/00b-add_author-{file.run_id}.log"
    shell:
        """
        python ./scripts/add_authors.py \
            {input.metadata} \
            {output.metadata} 2> {log}
        """

rule seq_index:
    message: "Create index for sequences"
    input:
        sequences = rules.parse.output.sequences
    output:
        index = "results/sequences.fasta.idx"
    log:
        "logs/01a-seq_index-{file.run_id}.log"
    shell:
        """
        augur index --sequences {input.sequences} \
            --output {output.index} 2> {log}
        """

rule filter_data:
    message: "Filter sequences by inclusion / exclusion criteria for {input.metadata} -> {output.metadata} "
    input:
        sequences = rules.parse.output.sequences,
	metadata = rules.add_authors.output.metadata,
        index = rules.seq_index.output.index,
        include = file.include,
        exclude = file.exclude
    output:
        sequences = "results/sequences.fasta",
	metadata = "results/metadata.tsv"
    log:
        "logs/01b-filter_data-{file.run_id}.log"
    shell:
        """
        augur filter --metadata {input.metadata} \
            --sequences {input.sequences} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --exclude {input.exclude} \
            --include {input.include} 2> {log}
        """

rule create_colors:
    message:
        "Creating custom color scale in {output.colors}"
    input:
        metadata = rules.filter_data.output.metadata
    output:
        colors = "results/colors.tsv"
    log:
        "logs/02a-create_color-{file.run_id}.log"
    shell:
        """
        python ./scripts/make_colors.py {input.metadata} {output.colors} 2> {log}
        """

rule create_lat_longs:
    message:
        "Creating lat/longs in {output.lat_longs}"
    input:
        metadata = rules.filter_data.output.metadata
    output:
        lat_longs = "results/lat_longs.tsv"
    log:
        "logs/02b-create_lat_longs-{file.run_id}.log"
    shell:
        """
        python ./scripts/create_lat_longs.py {input.metadata} {output.lat_longs} 2> {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter_data.output.sequences,
        reference = file.reference
    output:
        alignment = "results/aligned.fasta"
    log:
        "logs/03a-align-{file.run_id}.log"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads 40 \
            --fill-gaps 2> {log}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    log:
        "logs/04a-tree-{file.run_id}.log"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method raxml \
#            --tree-builder-args '-m GTRCAT -p $RANDOM -f d -x $RANDOM -N 100' \
#            --override-default-args \
            --nthreads auto 2> {log}
        """

# add bootstrap x100 values on tree_raw, root at AY765264
# iqtree2 -s results/aligned.fasta --seqtype DNA -t results/tree_raw.nwk -B 1000 --prefix results/tree_bsv -T 20 --keep-ident -alninfo -o AY765264
rule vcf:
    message: "Create VCF files for MATree"
    input:
        alignment = rules.align.output.alignment
    output:
        vcf_file = "results/merged.vcf"
    params:
        reference = "NC_009942.1"
    shell:
        """
        faToVcf -ref={params.reference} \
            {input.alignment} \
            {output.vcf_file}
        """
rule matree:
    message:
        """
        Creating MATree
        """
    input:
        tree = rules.tree.output.tree,
        vcf_file = rules.vcf.output.vcf_file
    output:
        matree = "results/MATree.pb"
    shell:
        """
        usher -t {input.tree} \
            -v {input.vcf_file} \
            -o {output.matree} 
        """

rule root_tree:
    message:
        """
        Creating root_treeed tree
          - setting root of the tree as {params.root}
          - can use AF481864 for root, a pre-NY99 sequence as root for lineage 1a builds
          - OQ067500 is an outgroup virus - Koutango virus
        """
    input:
        matree = rules.matree.output.matree
    output:
        tree = "results/tree.nwk"
    params:
        root = "OQ067500.1"
    shell:
        """
        matUtils extract \
            -i {input.matree} \
            -y {params.root} \
            -t {output.tree} 
        """
		
rule annotate:
    message: "Annotate MATree with known lineage data"
    input:
        matree = rules.root_tree.output.matree,
        clades = file.known_clades
    output:
        matree = "MATree_annotated.pb",
        mutations = "MATree_annot_mutations.tsv",
        clade_nodes = "MATree_annot_clade-details.tsv"
    params:
        outdir = "./results/"
    shell:
        """
        matUtils annotate -i {input.matree} \
            -c {input.clades} \
            -o {output.matree} \
            -d {params.outdir} \
            -u {output.mutations} \
            -D {output.clade_nodes}
        """

## melt clade-details.tsv

rule lineage_base:
    message: "Retrieve annotated lineages"
    input:
        matree = rules.annotate.output.matree
    output:
        scores = "MATree_sample-scores.tsv",
		clades = "MATree_sample-clades.tsv",
        aberrant = "MATree_sample-aberrant.tsv"
    params:
        outdir = "./results/"
    shell:
        """
        matUtils summary -i {input.matree} \
            -d {params.outdir} \
            -s {output.scores} \
            -a {output.aberrant} \
            -C {output.clades}
        """

rule lineage_new:
    message: "Create new lineages based on existing lineages"
    input:
        matree = rules.annotate.output.matree
    output:
        matree = "results/MATree_lineages.pb",
        clades = "results/autolin_clades.tsv",
        labels = "results/autolin_labels.tsv"
    shell:
        """
        python3 ./scripts/propose_sublineages.py \
            -i {input.matree} \
            -r \
            -o {output.matree} \
            -d {output.clades} \
            -l {output.labels}
        """
		
rule annotate_lineage:
    message: "Annotate MATree with known lineage data"
    input:
        matree = rules.root_tree.output.matree,
        clades = rules.lineage_new.output.labels
    output:
        matree = "MATree_annotated.pb",
        mutations = "MATree_annot_mutations.tsv",
        clade_nodes = "MATree_annot_clade-details.tsv"
    params:
        outdir = "./results/"
    shell:
        """
        matUtils annotate -i {input.matree} \
            -c {input.clades} \
            -o {output.matree} \
            -d {params.outdir} \
			-f 0.95 \
            -u {output.mutations} \
            -D {output.clade_nodes}
        """
# UShER re-root, nwk re-root AY765264, NY99 `matUtil extract `

rule ancestral:
    message: 
        """
        Reconstructing ancestral sequences and mutations
          - using {input.root} for mutation calling
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.matree.output.tree,
        alignment = rules.align.output.alignment,
        root = file.reference
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    log:
        "logs/05a-ancestral-{file.run_id}.log"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --root-sequence {input.root} \ 
            --inference {params.inference} 2> {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.root_tree.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = file.reference
    output:
        node_data = "results/aa_muts.json"
    log:
        "logs/05b-translate-{file.run_id}.log"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}  2> {log}
        """

rule clades:
    message: "Setting clade membership using clade defining mutations"
    input:
        tree = rules.root_tree.output.tree,
        aa_nodes = rules.translate.output.node_data,
        aa_clades = file.clades
    output:
        node_data = "results/clade_membership.json",
    log:
        "logs/06a-clades-{file.run_id}.log"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.aa_nodes} \
            --clades {input.aa_clades} \
            --output {output.node_data} 2> {log}
        """

rule usher_clades:
    message: "Setting strain membership using clade defining mutations"
    input:
        tree = rules.root_tree.output.tree,
        nt_nodes = rules.ancestral.output.node_data,
        nt_strains = file.usher_mutations
    output:
        node_data = "results/usher_clades.json"
    params:
        labels = "usher_lineages"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_nodes} \
            --clades {input.nt_strains} \
            --output {output.node_data} \
            --membership-name {params.labels} 
        """


rule lineage_clades:
    message: "Setting strain membership using clade defining mutations"
    input:
        tree = rules.root_tree.output.tree,
        nt_nodes = rules.ancestral.output.node_data,
        nt_strains = file.lineage_mutations
    output:
        node_data = "results/lineage_clades.json"
    params:
        labels = "autolin_lineages"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_nodes} \
            --clades {input.nt_strains} \
            --output {output.node_data} \
            --membership-name {params.labels} 
        """
#!!!rule autolin-clades:
#autolin.pb clade.tsv (sample-newclades)
#matUtils annotate final.pb - autolin_details.tsv --> melt autolin_mutations 
# nextstrain clades --> autolin_lineages.json

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.root_tree.output.tree,
        metadata = rules.filter_data.output.metadata
    output:
        node_data = "results/traits.json"
    params:
        columns = "country clade_membership strains lineages"
    log:
        "logs/07a-traits-{file.run_id}.log"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence 2> {log}
        """


## merge metadata script -> "results/metadata_fin.tsv"

rule export_nextclade:
    message:
        """
        Exporting data files for for auspice using V2 JSON schema with all lineages and inferred lineages
        Including
          - nucleotide {input.nt_muts}
          - amino acid {input.aa_muts}
          - using {rules.ancestral.params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.root_tree.output.tree,
        metadata = "results/metadata_fin.tsv",
        clades = rules.clades.output.node_data,
        strains = rules.usher_clades.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = rules.create_colors.output.colors,
        lat_longs = rules.create_lat_longs.output.lat_longs,
        auspice_config = "config/auspice_config_v2.json"
    output:
        auspice = "auspice/tree.json"
    log:
        "logs/10c-export_nextclade-{all.output.run_id}.log"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.nt_muts} {input.aa_muts} {input.clades} {input.traits} {input.strains} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice} 2> {log}
        """


rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
