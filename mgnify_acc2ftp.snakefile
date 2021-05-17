import pandas as pd

#m = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv", sep = "\t")
#m = m.sample(n = 1000)
m = pd.read_csv("tmp-1000r-mgnify-genomes-all_metadata.tsv", sep = "\t")
m = m.head(n = 100)
m.dropna(subset = ["Sample_accession"], inplace = True)
GENOMES = m['Genome'].unique().tolist()
ACCESSIONS = m['Sample_accession'].unique().tolist()

rule all:
    input:
        expand("inputs/mgnify_raw_reads_links/{accession}.tsv", accession = ACCESSIONS)

rule refinem_download_ftp_links_for_raw_fastq:
    output: "inputs/mgnify_raw_reads_links/{accession}.tsv"
    params: url = lambda w: "'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + w.accession + "&result=read_run&fields=fastq_ftp'"
    shell:'''
    wget -O {output} {params.url}
    '''

