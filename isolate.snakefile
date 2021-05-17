import pandas as pd
import re

#gtdb_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_taxonomy_r89.tsv"
gtdb_url = "inputs/bac120_taxonomy_r89.tsv"
gtdb = pd.read_csv(gtdb_url, sep = "\t")
GENOMES = gtdb[gtdb.columns[0]].unique().tolist()
GENOMES = [genome for genome in GENOMES if 'RS_' in genome] # filter to refseq
GENOMES = [re.sub("RS_", "", genome) for genome in GENOMES]

gtdb_failed = pd.read_csv("failed_gtdb_140k_download.txt", header = None)
gtdb_failed = gtdb_failed[gtdb_failed.columns[0]].unique().tolist()
GENOMES = set(GENOMES) - set(gtdb_failed) 

rule all:
    input: 
        #expand("inputs/refseq_gtdb_isolates/{genome}.fna.gz", genome = GENOMES)
        expand("outputs/charcoal_isolate/{genome}.fa.clean.fa.gz", genome = GENOMES),

rule download_datasets_tool:
    output: "scripts/datasets"
    shell:'''
    # specific to macs
    # wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
    chmod +x {output}
    '''

rule datasets_download:
    input: "scripts/datasets"
    output: "inputs/refseq_gtdb_isolates_zip/{genome}.zip"
    params: genome = lambda wildcards: wildcards.genome
    shell:'''
    scripts/datasets download assembly {params.genome} -f {output}   
    '''

rule unzip_genomes:
    input: "inputs/refseq_gtdb_isolates_zip/{genome}.zip"
    output: "inputs/refseq_gtdb_isolates/{genome}.fna.gz"
    shell:'''
    unzip -p {input} ncbi_dataset/data/{wildcards.genome}/*.fna | gzip > {output}
    '''
#rule download_genomes:
#    output: "inputs/refseq_gtdb_isolates/{genome}_genomic_refseq.fna.gz"
#    params: genome = lambda wildcards: wildcards.genome
#    #conda: "envs/biomartr.yml"
#    script: "scripts/download_gtdb_genomes_biomartr.R"

rule make_mag_list: 
    """
    This rule will over specify the number of files that charcoal needs to run on if there are extra genome files in this directory.
    """
    input: expand("inputs/refseq_gtdb_isolates/{genome}.fna.gz", genome = GENOMES)
    output: "inputs/charcoal_conf/isolate_genome_list.txt"
    params: indir = "inputs/refseq_gtdb_isolates"
    shell:'''
    ls {params.indir} > {output}
    '''

rule run_charcoal:
    input: 
        genomes = expand("inputs/refseq_gtdb_isolates/{genome}.fna.gz", genome = GENOMES),
        mag_list = "inputs/charcoal_conf/isolate_genome_list.txt",
        conf = "inputs/charcoal_conf/isolate_genomes.conf",
    output: 
        expand("outputs/charcoal_isolate/{genome}.fa.clean.fa.gz", genome = GENOMES),
        expand("outputs/charcoal_isolate/{genome}.fa.dirty.fa.gz", genome = GENOMES)
    conda: "envs/charcoal.yml"
    benchmark: "benchmarks/charcoal.benchmark"
    shell:'''
    charcoal run {input.conf} -j 16 --nolock --no-use-conda
    '''
