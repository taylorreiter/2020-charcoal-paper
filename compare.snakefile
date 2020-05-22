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
        expand("outputs/charcoal/{genome}.fa.clean.fa", genome = GENOMES)
        

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
    gunzip {input}
    '''

rule run_checkm_lineage_qf_clean:
    """
    Reports general summary statistics (% complete, % contamination) for genome (bins).
    """
    input: expand("outputs/charcoal/{genome}.fa.clean.fa.gz", genome = GENOMES)
    output: "outputs/charcoal_clean_checkm/completeness.tsv"
    params: 
        indir = "outputs/charcoal",
        threads = "4",
        outdir = "outputs/charcoal_clean_checkm/"
    conda: "envs/checkm.yml"
    #benchmark: "benchmarks/{genome}_checkm_clean.benchmark"
    shell:'''
    checkm lineage_wf \
        --file {output} \
        --tab_table \
        --extension .fa.clean.fa.gz \
        --threads {params.threads} \
        {params.indir} {params.outdir} 
    '''     

rule run_checkm_qa_clean:
    """
    Reports name of contig on which contaminant marker gene resides, as well as marker gene identity.
    """
    input: marker_file = "charcoal_clean_checkm/lineage.ms"
    output: ""
    conda: "envs/checkm.yml"
    benchmark: "benchmarks/checkm_qa_clean.benchmark"
    params: 
        threads = "4",
        indir = "charcoal_clean_checkm" 
    shell:'''
    checkm qa -o 8 -f {output} --tab_table -t {params.threads} {input.marker_file} {params.indir}
    '''

rule run_prokka_dirty:
    input: "outputs/charcoal/{genome}.fa.dirty.fa.gz"
    output: "outputs/charcoal_dirty_prokka/{genome}.fna"
    params: outdir = "outputs/charcoal_dirty_prokka" 
    conda: "envs/prokka.yml"
    benchmark: "benchmarks/{genome}_prokka_dirty.benchmark"
    shell:'''
    prokka {input} --outdir {params.outdir} --prefix {wildcards.genome} --metagenome --force --locustag {wildcards.genome}
    '''

############################################
### Compare against MAGpurify
############################################


rule run_magpurify_phylo_markers:
    shell:'''
    magpurify phylo-markers example/test.fna example/output
    '''

rule run_magpurify_clade_markers:
    shell:'''
    magpurify clade-markers example/test.fna example/output
    '''

rule run_magpurify_tetra_freq:
    shell:'''
    magpurify tetra-freq example/test.fna example/output
    '''

rule run_magpurify_gc_content:
    shell:'''
    magpurify gc-content example/test.fna example/output
    '''

rule run_magpurify_known_contam:
    shell:'''
    magpurify known-contam example/test.fna example/output
    '''

rule run_magpurify_clean_bin:
    shell:'''
    magpurify clean-bin example/test.fna example/output example/test_cleaned.fna
    '''

#############################################
### Compare against refineM
#############################################

rule run_refineM:
