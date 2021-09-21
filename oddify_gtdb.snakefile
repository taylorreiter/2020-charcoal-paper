# download gtdb rs202 lca db
# 

DATABASES = ['gtdb-rs202.genomic', 'gtdb-rs202.genomic-reps']
KSIZE = [31, 51]

rule all:
    input:
        expand("outputs/{db}-oddities-k{ksize}.csv", db = DATABASES, ksize = KSIZE)
         
rule download_lca_dbs:
    output: expand("inputs/lca_dbs/{db}.k{ksize}.lca.json.gz", db = DATABASES, ksize = KSIZE)
    shell:'''
    wget -O inputs/lca_dbs/gtdb-rs202.genomic.k31.lca.json.gz https://osf.io/9xdg2/download
    wget -O inputs/lca_dbs/gtdb-rs202.genomic.k51.lca.json.gz https://osf.io/3cdp6/download
    wget -O inputs/lca_dbs/gtdb-rs202.genomic-reps.k31.lca.json.gz https://osf.io/ypsjq/download
    wget -O inputs/lca_dbs/gtdb-rs202.genomic-reps.k51.lca.json.gz https://osf.io/297dp/download
    '''

rule download_oddify_scripts:
    output:
        find = "scripts/find-oddities.py",
        examine = "scripts/find-oddities-examine.py"
    shell:'''
    wget -O {output.find} https://raw.githubusercontent.com/dib-lab/sourmash-oddify/eab9a2118246c3ed2189fe166c6dc388a888ecee/scripts/find-oddities.py
    wget -O {output.examine} https://raw.githubusercontent.com/dib-lab/sourmash-oddify/update_lca_db_use/scripts/find-oddities-examine.py
    '''

rule find_oddities_txt:
    input:
        lca="inputs/lca_dbs/{db}.k{ksize}.lca.json.gz",
        find = "scripts/find-oddities.py"
    output:
        "outputs/{db}-oddities-k{ksize}.txt",
        "outputs/{db}-oddities-k{ksize}.csv"
    params:
        outdir="outputs"
    conda: "envs/sourmash.yml"
    resources: mem_mb = 32000
    shell:
        "scripts/find-oddities.py {input.lca} --lowest-rank=superkingdom --minimum-hashes=10 --prefix={params.outdir}/{wildcards.db}-oddities-k{wildcards.ksize} > {output[0]}"

rule make_oddities_examine_txt:
    input:
        "outputs/{db}-oddities-k{ksize}.txt",
        "outputs/{db}-genomes-k{ksize}.lca.json.gz",
        "outputs/{db}-oddities-k{ksize}.csv",
        "scripts/find-oddities-examine.py"
    output:
        "outputs/{db}-oddities-k{ksize}.examine.txt"
    params:
        outdir="outputs",
        genomes_dir="/home/ntpierce/2021-rank-compare/genbank/genomes",
        extension=".fna.gz"
    resources: mem_mb = 8000
    conda: "envs/sourmash.yml"
    shell:
        "scripts/find-oddities-examine.py {params.outdir}/{wildcards.db}-oddities-k{wildcards.ksize}.csv {params.genomes_dir} --percent-threshold=95 --length-threshold=0 --genome-extension={params.extension} > {output}"
