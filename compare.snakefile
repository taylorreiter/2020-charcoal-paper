import pandas as pd
import random 

#m = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv", sep = "\t")
#m = m.sample(n = 1000)
m = pd.read_csv("tmp-1000r-mgnify-genomes-all_metadata.tsv", sep = "\t")
m = m.head(n = 100)
GENOMES = m['Genome'].unique().tolist()
METAGENOMES = m['Sample_accession'].unique().tolist()

rule all:
    input:
        #expand("inputs/mgnify_genomes/human-gut/v1.0/{genome}.gff3.gz", genome = GENOMES)
        expand("outputs/charcoal/{genome}.fa.clean.fa", genome = GENOMES),
        "outputs/charcoal_clean_checkm_qa/qa.tsv",
        expand("outputs/charcoal_dirty_prokka/{genome}.fna", genome = GENOMES),
        "outputs/charcoal_dirty_checkm_qa/qa.tsv",
        expand("outputs/magpurify_clean/{genome}_magpurify.fna", genome = GENOMES)
        

rule download_genomes: 
    output: 
        genome="inputs/mgnify_genomes/human-gut/v1.0/{genome}.gff3.gz",
    run:
        row = m.loc[m['Genome'] == wildcards.genome]
        gff = row['FTP_download'].values
        gff = gff[0]
        shell("wget -O {output.genome} {gff}")


rule unembed_gff_fasta:
    input: "inputs/mgnify_genomes/human-gut/v1.0/{genome}.gff3.gz"
    output: "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa"
    shell:"""
    zcat {input} | awk 'BEGIN {{ doprint = 0}}; \
                        doprint == 1 {{ print $0 }}; \
                        $0 ~ /#FASTA/ {{ doprint = 1 }}' > {output}
    """

rule make_mag_list: 
    """
    This rule will over specify the number of files that charcoal needs to run on if there are extra genome files in this directory.
    """
    input: expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES)
    output: "inputs/charcoal_conf/mgnify_genome_list.txt"
    params: indir = "outputs/mgnify_genomes/human-gut/v1.0"
    shell:'''
    ls {params.indir} > {output}
    '''

rule run_charcoal:
    input:
        genomes = expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
        mag_list = "inputs/charcoal_conf/mgnify_genome_list.txt",
        conf = "inputs/charcoal_conf/mgnify_genomes.conf",
    output: 
        expand("outputs/charcoal/{genome}.fa.clean.fa.gz", genome = GENOMES),
        expand("outputs/charcoal/{genome}.fa.dirty.fa.gz", genome = GENOMES)
    conda: "envs/charcoal.yml"
    benchmark: "benchmarks/charcoal.benchmark"
    shell:'''
    charcoal run {input.conf} -j 4 --nolock --no-use-conda
    '''

rule gunzip_charcoal_clean:
    input: "outputs/charcoal/{genome}.fa.clean.fa.gz"
    output: "outputs/charcoal/{genome}.fa.clean.fa"
    shell:'''
    if [ -s {input} ]; then
        gunzip {input}
    else
        mv {input} {output}
    fi
    '''

rule run_checkm_lineage_qf_clean:
    """
    Reports general summary statistics (% complete, % contamination) for genome (bins).
    """
    input: expand("outputs/charcoal/{genome}.fa.clean.fa", genome = GENOMES)
    output: 
        comp = "outputs/charcoal_clean_checkm/completeness.tsv",
        lin = "outputs/charcoal_clean_checkm/lineage.ms"
    params: 
        indir = "outputs/charcoal",
        threads = "4",
        outdir = "outputs/charcoal_clean_checkm/"
    conda: "envs/checkm.yml"
    #benchmark: "benchmarks/{genome}_checkm_clean.benchmark"
    shell:'''
    checkm lineage_wf \
        --file {output.comp} \
        --tab_table \
        -x .fa.clean.fa \
        --threads {params.threads} \
        {params.indir} {params.outdir} 
    '''     

rule run_checkm_qa_clean:
    """
    Reports name of contig on which contaminant marker gene resides, as well as marker gene identity.
    """
    input: marker_file = "outputs/charcoal_clean_checkm/lineage.ms"
    output: "outputs/charcoal_clean_checkm_qa/qa.tsv"
    conda: "envs/checkm.yml"
    benchmark: "benchmarks/checkm_qa_clean.benchmark"
    params: 
        threads = "4",
        indir = "outputs/charcoal_clean_checkm" 
    shell:'''
    checkm qa -o 8 -f {output} --tab_table -t {params.threads} {input.marker_file} {params.indir}
    '''

rule gunzip_charcoal_dirty:
    input: "outputs/charcoal/{genome}.fa.dirty.fa.gz"
    output: "outputs/charcoal/{genome}.fa.dirty.fa"
    shell:'''
    if [ -s {input} ]; then
        gunzip {input}
    else
        mv {input} {output}
    fi
    '''

rule run_prokka_dirty:
    input: "outputs/charcoal/{genome}.fa.dirty.fa"
    output: "outputs/charcoal_dirty_prokka/{genome}.fna"
    params: outdir = "outputs/charcoal_dirty_prokka" 
    conda: "envs/prokka.yml"
    benchmark: "benchmarks/{genome}_prokka_dirty.benchmark"
    shell:'''
    if [ -s {input} ]; then
        prokka {input} --outdir {params.outdir} --prefix {wildcards.genome} --metagenome --force --locustag {wildcards.genome}
    else
        touch {output}
    fi
    '''

rule run_checkm_lineage_qf_dirty:
    """
    Reports general summary statistics (% complete, % contamination) for genome (bins).
    """
    input: expand("outputs/charcoal/{genome}.fa.dirty.fa", genome = GENOMES)
    output: 
        comp = "outputs/charcoal_dirty_checkm/completeness.tsv",
        lin = "outputs/charcoal_dirty_checkm/lineage.ms"
    params: 
        indir = "outputs/charcoal",
        threads = "4",
        outdir = "outputs/charcoal_dirty_checkm/"
    conda: "envs/checkm.yml"
    #benchmark: "benchmarks/{genome}_checkm_clean.benchmark"
    shell:'''
    checkm lineage_wf \
        --file {output.comp} \
        --tab_table \
        -x .fa.dirty.fa \
        --threads {params.threads} \
        {params.indir} {params.outdir} 
    '''     

rule run_checkm_qa_dirty:
    """
    Reports name of contig on which contaminant marker gene resides, as well as marker gene identity.
    """
    input: marker_file = "outputs/charcoal_dirty_checkm/lineage.ms"
    output: "outputs/charcoal_dirty_checkm_qa/qa.tsv"
    conda: "envs/checkm.yml"
    benchmark: "benchmarks/checkm_qa_dirty.benchmark"
    params: 
        threads = "4",
        indir = "outputs/charcoal_dirty_checkm" 
    shell:'''
    checkm qa -o 8 -f {output} --tab_table -t {params.threads} {input.marker_file} {params.indir}
    '''

############################################
### Compare against MAGpurify
############################################

rule download_magpurify_db:
    output: "inputs/magpurify_db/MAGpurify-db-v1.0.tar.bz2"
    shell:'''
    wget -O {output} https://zenodo.org/record/3688811/files/MAGpurify-db-v1.0.tar.bz2
    '''

rule decompress_magpurify_db:
    input: "inputs/magpurify_db/MAGpurify-db-v1.0.tar.bz2"
    output: "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm" 
    params: outdir = "inputs/magpurify_db"
    shell:'''
    tar -jxvf {input} -C {params.outdir}    
    '''

rule run_magpurify_phylo_markers:
    input: 
        genome =  "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db = "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm" 
    output: "outputs/magpurify/{genome}/phylo-markers/genes.out" 
    conda: "envs/magpurify.yml"
    params: outdir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    benchmark: "benchmarks/magpurify_phylo_markers_{genome}.benchmark"
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify phylo-markers {input.genome} {params.outdir}
    '''

rule run_magpurify_clade_markers:
    input: 
        genome = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db =  "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm" 
    output: "outputs/magpurify/{genome}/clade-markers/genes.out" 
    conda: "envs/magpurify.yml"
    benchmark: "benchmarks/magpurify_clade_markers_{genome}.benchmark"
    params: outdir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify clade-markers {input.genome} {params.outdir}
    '''

rule run_magpurify_tetra_freq:
    input: 
        genome = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db = "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm" 
    output: "outputs/magpurify/{genome}/tetra-freq/flagged_contigs" 
    params: outdir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    conda: "envs/magpurify.yml"
    benchmark: "benchmarks/magpurify_tetra_freq_{genome}.benchmark"
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify tetra-freq {input.genome} {params.outdir}
    '''

rule run_magpurify_gc_content:
    input: 
        genome = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db = "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm" 
    output: "outputs/magpurify/{genome}/gc-content/flagged_contigs" 
    params: outdir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    conda: "envs/magpurify.yml"
    benchmark: "benchmarks/magpurify_gc_content_{genome}.benchmark"
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify gc-content {input.genome} {params.outdir}
    '''

rule run_magpurify_known_contam:
    input: 
        genome = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db ="inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm",
    output: "outputs/magpurify/{genome}/known-contam/flagged_contigs", 
    params: outdir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    conda: "envs/magpurify.yml"
    benchmark: "benchmarks/magpurify_known_contam_{genome}.benchmark"
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify known-contam {input.genome} {params.outdir}
    '''

rule run_magpurify_clean_bin:
    input: 
        genome = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        db = "inputs/magpurify_db/MAGpurify-db-v1.0/phylo-markers/PhyEco.hmm", 
        gc = "outputs/magpurify/{genome}/gc-content/flagged_contigs", 
        tf = "outputs/magpurify/{genome}/tetra-freq/flagged_contigs", 
        cm = "outputs/magpurify/{genome}/clade-markers/genes.out", 
        pm = "outputs/magpurify/{genome}/phylo-markers/genes.out",
        kc = "outputs/magpurify/{genome}/known-contam/flagged_contigs", 
    output: "outputs/magpurify_clean/{genome}_magpurify.fna"
    conda: "envs/magpurify.yml"
    params: indir = lambda wildcards: "outputs/magpurify/" + wildcards.genome
    benchmark: "benchmarks/magpurify_clean_{genome}.benchmark"
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify clean-bin {input.genome} {params.indir} {output}
    '''

#############################################
### Compare against refineM
#############################################

rule refinem_download_db:
    output:""
    shell:'''
    wget -O {output}
    '''

rule refinem_download_raw_reads_forward:
    output: 
        r1="inputs/raw_mgnify_reads/{metagenome}_1.fastq.gz",
    run:
        metagenome = wildcards.metagenome
        metagenome_root = metagenome[0:6]
        download_ftp = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/" + metagenome_root + metagenome + "_1.fastq.gz"
        shell("wget -O {output.r1} {fastq_1}")


rule refinem_download_raw_reads_reverse:
    output: 
        r2="inputs/raw_mgnify_reads/{metagenome}_2.fastq.gz",
    run:
        metagenome = wildcards.metagenome
        metagenome_root = metagenome[0:6]
        download_ftp = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/" + metagenome_root + metagenome + "_2.fastq.gz"
        shell("wget -O {output.r2} {fastq_2}")

rule refinem_index_genomes:
    input: 
        genome ="outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
    output: "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa.bwt"
    conda: "envs/refinem.yml"
    shell:'''
    bwa index {input}
    '''

rule refinem_align_reads:
    input:
        genome ="outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
        indx = "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa.bwt"
        r1 = expand("inputs/raw_mgnify_reads/{metagenome}_1.fastq.gz", metagenome = METAGENOMES),
        r2 = expand("inputs/raw_mgnify_reads/{metagenome}_2.fastq.gz", metagenome = METAGENOMES),
    output: "outputs/refinem/bams/{genome}.bam"
    benchmark: "benchmarks/refinem_{genome}_bwa_align_reads.txt"
    run:
        row = m.loc[m['Genome'] == wildcards.genome]
        metagenome = row['Sample_accession'].values
        metagenome = metagenome[0]
        r1 = "inputs/raw_mgnify_reads/" + metagenome + "_1.fastq.gz"
        r2 = "inputs/raw_mgnify_reads/" + metagenome + "_2.fastq.gz"
        shell("bwa mem {input.genome} {r1} {r2} | samtools sort -o {output} -")


rule refinem_scaffold_stats:
    input: 
        expand("outputs/refinem/bams/{genome}.bam", genome = GENOMES),
        expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
    output: "outputs/refinem/scaffold_stats/scaffold_stats.tsv" 
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        bamdir = "outputs/refinem/bams", 
        outdir = "outputs/refinem/scaffold_stats"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_scaffold_stats.txt"
    shell:'''
    refinem scaffold_stats -c 16 scaffold_stats.tsv {params.bindir} {params.outdir} {params.bamdir}
    '''

rule refinem_outliers:
    input: "outputs/refinem/scaffold_stats/scaffold_stats.tsv" 
    output: "outputs/refinem/outliers/outliers.tsv"
    params: outdir = "outputs/refinem/outliers"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_outliers.txt"
    shell:'''
    refinem outliers {input} <outlier_output_dir>
    '''

rule refinem_filter_bins_stats:
    input: "outputs/refinem/outliers/outliers.tsv"
    output: "outputs/refinem/filter_bin_stats/" # need to decide if I want all bins here or just modified bins. if mod, probs directory().
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        outdir = "outputs/refinem/filter_bins_stats/"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_filter_bins_stats.txt"
    shell:'''
    refinem filter_bins {params.bindir} {input} {params.outdir} # --modified_only will make this rule only output modified bins
    '''

rule refinem_call_genes:
    input: 
        expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
    output: # NEED TO FILL IN OUTPUT FILE (run first, to determine name of an output file)
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        genedir = "outputs/refinem/call_genes/"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_call_genes.txt"
    shell:'''
    refinem call_genes -c 40 {params.bindir} {params.genedir}
    '''

rule refinem_taxon_profile:
    input: 
        scaffold_stats = "outputs/refinem/scaffold_stats/scaffold_stats.tsv",
        db = ,# NEED TO FILL IN WITH DB DOWNLOAD
        tax = # NEED TO FILL IN WITH REF TAX DOWNLOAD
    output: "outputs/refinem/taxon_profile/" # NEED TO ADD A SPECIFIC FILE NAME
    params:
        genedir = "outputs/refinem/call_genes",
        outdir = "outputs/refinem/taxon_profile"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_taxon_profile.txt"
    shell:'''
    refinem taxon_profile -c 40 {params.genedir} {input.scaffold_stats} {input.db} {input.tax} {params.outdir}
    '''

rule refinem_taxon_filter:
    input: "outputs/refinem/taxon_profile/" # NEED TO ADD A SPECIFIC FILE NAME TO MATCH OUTPUT OF ABOVE RULE 
    output: "outputs/refinem/taxon_profile/taxon_filter.tsv" 
    params: indir = "outputs/refinem/taxon_profile"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_taxon_filter.txt"
    shell:'''
    refinem taxon_filter -c 40 {params.indir} {output}
    '''

rule refinem_filter_bins_taxon:
    input: 
        tax_filt = "outputs/refinem/taxon_profile/taxon_filter.tsv",
        genome = expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
    output: "outputs/refinem/filter_bins_taxon" # NEED TO ADD SPECIFIC FILE OUTPUTS
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        outdir = "outputs/refinem/filter_bins/taxon"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_filter_bins_taxon.txt"
    shell: '''
    refinem filter_bins {params.bindir} {input.tax_filt} {params.outdir}
    '''
