import csv
import os
import pandas as pd

tax_info = pd.read_csv("inputs/gtdb-rs202.taxonomy.with-repinfo.csv", header=0)
tax_info['genome_fastafile'] = 'genbank/genomes/'+ tax_info['ident'] + "_genomic.fna.gz" # add intended fastapaths

ACCS =  tax_info["ident"]
#tax_info.set_index("ident", inplace=True)

rule all:
    input:
        expand("outputs/gtdb_rs202_genomes/genomes/{acc}_genomic.fna.gz", acc=ACCS)

# ad genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'outputs/gtdb_rs202_genomes/info/{acc}.info.csv'
    shell: """
        python -Werror genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """

# download actual genomes!
rule download_matching_genome_wc:
     input:
         csvfile = ancient('outputs/gtdb_rs202_genomes/info/{acc}.info.csv')
     output:
         genome = "outputs/gtdb_rs202_genomes/genomes/{acc}_genomic.fna.gz"
     run:
         with open(input.csvfile, 'rt') as infp:
             r = csv.DictReader(infp)
             rows = list(r)
             assert len(rows) == 1
             row = rows[0]
             acc = row['acc']
             assert wildcards.acc.startswith(acc)
             url = row['genome_url']
             name = row['ncbi_tax_name']

             print(f"downloading genome for acc {acc}/{name} from NCBI...",
                   file=sys.stderr)
             with open(output.genome, 'wb') as outfp:
                 with urllib.request.urlopen(url) as response:
                     content = response.read()
                     outfp.write(content)
                     print(f"...wrote {len(content)} bytes to {output.genome}",
                           file=sys.stderr)

rule make_charcoal_genome_list: 
    """
    This rule will over specify the number of files that charcoal needs to run on if there are extra genome files in this directory.
    """
    input: expand("outputs/gtdb_rs202_genomes/genomes/{acc}_genomic.fna.gz", acc = ACC)
    output: "outputs/charcoal_conf/gtdb_rs202_genome_list.txt"
    params: indir = "outputs/gtdb_rs202_genomes/genomes"
    shell:'''
    ls {params.indir}/*gz | xargs -n 1 basename > {output}
    '''

rule make_charcoal_lineage_sheet:

rule run_charcoal:
    input: 
        genomes = expand("outputs/gtdb_rs202_genomes/genomes/{acc}_genomic.fna.gz", acc = ACC),
        genome_list = "outputs/charcoal_conf/gtdb_rs202_genome_list.txt",
        conf = "inputs/charcoal_conf/gtdb_rs202_genomes.conf",
    output: 
        expand("outputs/gtdb_rs202_charcoal/{genome}.fa.clean.fa.gz", acc = ACC),
        expand("outputs/gtdb_rs202_charcoal/{genome}.fa.dirty.fa.gz", acc = ACC)
    conda: "envs/charcoal.yml"
    benchmark: "benchmarks/charcoal_gtdb_rs202.benchmark"
    shell:'''
    charcoal run {input.conf} -j 16 --nolock --no-use-conda
    '''
