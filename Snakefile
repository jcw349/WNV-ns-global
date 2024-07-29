rule all:
    input:
        auspice_tree = "auspice/WNV-global.json"
        auspice_tree_inferred = "auspice/WNV-global-infer.json"
        auspice_tree_1 = "auspice/WNV-1.json"
	auspice_tree_2 = "auspice/WNV-2.json"

rule files:
    params:
        input_fasta = "data/full_dataset.fasta",
        input_metadata = "data/headers.csv",
        reference = "config/reference.gb",
        auspice_config = "config/auspice_config_v2.json",
        clades = "config/clade.tsv",
        strains = "config/strain.tsv
        lat_longs = "config/lat_longs.tsv",
        exclude= "config/exclude.txt",
        include= "config/include.txt"

files = rules.files.params

rule parse:
    message:
        "Parsing {input.sequences}, {input.metadata} and forming FASTA + metadata TSV"
    input:
        sequences = files.input_fasta,
        metadata = files.input_metadata
    output:
        sequences = "results/sequences_raw.fasta",
        metadata = "results/metadata_sans_authors.tsv"
    shell:
        """
        python ./scripts/parse_fasta_csv.py \
            {input.sequences} \
            {input.metadata} \
            {output.sequences} \
            {output.metadata}
        """

rule add_authors:
    message:
        "Adding authors to {input.metadata} -> {output.metadata} by collecting info from ENTREZ"
    input:
        metadata = rules.parse.output.metadata
    output:
        metadata = "results/metadata_raw.tsv"
    shell:
        """
        python ./scripts/add_authors.py \
            {input.metadata} \
            {output.metadata}
        """

rule seq_index:
    message: "Create index for sequences"
    input:
        sequences = rules.parse.output.sequences
    output:
        index = "results/sequences.fasta.idx"
    shell:
        """
        augur index --sequences {input.sequences} \
            --output {output.index}
        """

rule filter_data:
    message: "Filter sequences by inclusion / exclusion criteria for {input.metadata} -> {output.metadata} "
    input:
        sequences = rules.parse.output.sequences
	metadata = rules.add_authors.output.metadata
        index = rule.seq_index.output.index
        include = files.include
        exclude = files.exclude
    output:
        sequences = "results/sequences.fasta"
	metadata = "results/metadata.tsv"
    shell:
        """
        augur filter --metadata {input.metadata} \
            --sequences {input.sequences} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --exclude {input.exclude} \
            --include {input.include} 
        """

rule create_colors:
    message:
        "Creating custom color scale in {output.colors}"
    input:
        metadata = rules.filter_data.output.metadata,
    output:
        colors = "results/colors.tsv"
    shell:
        """
        python ./scripts/make_colors.py {input.metadata} {output.colors}
        """

rule create_lat_longs:
    message:
        "Creating lat/longs in {output.lat_longs}"
    input:
        metadata = rules.filter_data.output.metadata,
    output:
        lat_longs = "results/lat_longs.tsv"
    shell:
        """
        python ./scripts/create_lat_longs.py {input.metadata} {output.lat_longs}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter_data.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method raxml \
            --nthreads auto
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
          - setting root of the tree as {params.root}, lineage 3 sequence--the most distantly related sequence from the rest of the tree
          - can use AF481864 for root, a pre-NY99 sequence as root for lineage 1a builds
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.filter_data.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        root = "AY765264" 
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --root {params.root}
        """
###
# ADD UShER and autolin
###
rule ancestral:
    message: 
        """
        Reconstructing ancestral sequences and mutations
          - using {input.root} for mutation calling
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
        root = rules.file.reference
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --root-sequence {input.root} \ 
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} 
        """
rule clades:
    message: "Setting clade membership using clade defining mutations"
    input:
        tree = rules.refine.output.tree,
        aa_nodes = rules.translate.output.node_data,
        aa_clades = files.clades
    output:
        node_data = "results/clade_membership.json",
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.aa_nodes} \
            --clades {input.aa_clades} \
            --output {output.node_data} 
        """

rule strains:
    message: "Setting strain membership using clade defining mutations"
    input:
        tree = rules.refine.output.tree,
        aa_nodes = rules.translate.output.node_data,
        aa_strains = files.strains
    output:
        node_data = "results/strain_membership.json",
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.aa_nodes} \
            --clades {input.aa_strains} \
            --output {output.node_data} 
        """
rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter_data.output.metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "country clade_membership strains lineages"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export_basic:
    message:
        """
        Exporting data files for for auspice using V2 JSON schema with all lineages
        Including
          - branch length {input.branch_lengths}
          - nucleotide {input.nt_muts}
          - amino acid {input.aa_muts}
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_authors.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = rules.create_colors.output.colors,
        lat_longs = rules.create_lat_longs.output.lat_longs,
        auspice_config = "config/auspice_config_v2.json"
    output:
        auspice = rule.all.auspice_tree
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice}
        """

rule export_inferred:
    message:
        """
        Exporting data files for for auspice using V2 JSON schema with all lineages and inferred lineages
        Including
          - branch length {input.branch_lengths}
          - nucleotide {input.nt_muts}
          - amino acid {input.aa_muts}
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.add_authors.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.node_data,
        strains = rules.strains.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = rules.create_colors.output.colors,
        lat_longs = rules.create_lat_longs.output.lat_longs,
        auspice_config = "config/auspice_config_v2.json"
    output:
        auspice = rule.all.auspice_tree_inferred
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} {input.strains} {input.traits} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice}
        """
rule export_I:
    message:
        """
        Exporting data files for for auspice using V2 JSON schema with lineage I and inferred lineage I
        Including
          - branch length {input.branch_lengths}
          - nucleotide {input.nt_muts}
          - amino acid {input.aa_muts}
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.refined_I.output.tree,
        metadata = rules.add_authors.output.metadata,
        branch_lengths = rules.refined_I.output.node_data,
        strains = rules.strains_I.output.node_data,
	lineages = rules.lineages_I.output.node_data,
        traits = rules.traits_I.output.node_data,
        nt_muts = rules.ancestral_I.output.node_data,
        aa_muts = rules.translate_I.output.node_data,
        colors = rules.create_colors.output.colors,
        lat_longs = rules.create_lat_longs.output.lat_longs,
        auspice_config = "config/auspice_config_v2.json"
    output:
        auspice = rule.all.auspice_tree_1
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} {input.strains} {input.traits} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice}
        """
rule export_II:
    message:
        """
        Exporting data files for for auspice using V2 JSON schema with lineage I and inferred lineage I
        Including
          - branch length {input.branch_lengths}
          - nucleotide {input.nt_muts}
          - amino acid {input.aa_muts}
          - using {params.inference} to infer ancestral maximum likelihood ancestral sequence states
        """
    input:
        tree = rules.refined_II.output.tree,
        metadata = rules.add_authors.output.metadata,
        branch_lengths = rules.refined_II.output.node_data,
        strains = rules.strains_II.output.node_data,
	lineages = rules.lineages_II.output.node_data,
        traits = rules.traits_II.output.node_data,
        nt_muts = rules.ancestral_II.output.node_data,
        aa_muts = rules.translate_II.output.node_data,
        colors = rules.create_colors.output.colors,
        lat_longs = rules.create_lat_longs.output.lat_longs,
        auspice_config = "config/auspice_config_v2.json"
    output:
        auspice = rule.all.auspice_tree_2
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} {input.strains} {input.traits} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --output {output.auspice}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
