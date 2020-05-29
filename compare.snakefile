import pandas as pd
import random 

#m = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv", sep = "\t")
#m = m.sample(n = 1000)
m = pd.read_csv("tmp-1000r-mgnify-genomes-all_metadata.tsv", sep = "\t")
m = m.head(n = 100)
GENOMES = m['Genome'].unique().tolist()

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
    shell:'''
    export MAGPURIFYDB=inputs/magpurify_db/MAGpurify-db-v1.0
    magpurify clean-bin {input.genome} {params.indir} {output}
    '''

#############################################
### Compare against refineM
#############################################

rule run_refineM:
