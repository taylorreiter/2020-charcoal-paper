import pandas as pd
import re

gtdb_url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_taxonomy_r89.tsv"
gtdb = pd.read_csv(gtdb_url, sep = "\t")
GENOMES = gtdb[gtdb.columns[0]].unique().tolist()
GENOMES = [genome for genome in GENOMES if 'RS_' in genome] # filter to refseq
GENOMES = [re.sub("RS_", "", genome) for genome in GENOMES]
GENOMES = GENOMES[0:5]

rule all:
    input: 
        expand("inputs/refseq_gtdb_isolates/{genome}.fna.gz", genome = GENOMES)

rule download_datasets_tool:
    output: "scripts/datasets"
    shell:'''
    # specific to macs
    wget -O {output} wget -O scripts/datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets
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

