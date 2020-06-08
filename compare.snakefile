import os
import pandas as pd
import random 

#m = pd.read_csv("http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv", sep = "\t")
#m = m.sample(n = 1000)
m = pd.read_csv("tmp-1000r-mgnify-genomes-all_metadata.tsv", sep = "\t", index_col=0)
m = m.head(n = 100)
GENOMES = m['Genome'].unique().tolist()
m.dropna(subset = ["Sample_accession"], inplace = True)
ACCESSIONS = m['Sample_accession'].unique().tolist()


##### refinem processing ######

## Establish infrastructure to link genome accession to
## metagenome read accessions. 

# read all into one pd dataframe, creating an "accession" column from the filename
# drop genomes that don't have an accession 
acc2ftp_files = expand("inputs/mgnify_raw_reads_links/{accession}.tsv", accession=ACCESSIONS)
acc2ftpDF = pd.concat([pd.read_csv(tsv, sep="\t").assign(Sample_accession=os.path.basename(tsv).rsplit(".tsv")[0]) for tsv in acc2ftp_files], ignore_index=True)
acc2ftpDF["ftp_links"] = acc2ftpDF["fastq_ftp"].str.split(";")

def generate_fastq_acc_and_output(row):
    only_or_second_filename = row["fastq_ftp"].rsplit("/", 1)[1]
    # get fastq_acc and libtype from the ftp filename
    if "_2" in only_or_second_filename:
    # pe: _1.fastq.gz, _2.fastq.gz ; SE: #.fastq.gz
        row["libtype"] = "paired"
        row["fastq_acc"] = only_or_second_filename.rsplit("_2.fastq.gz")[0]
        # generate output names from links, for use later
        row["outfile1"] = only_or_second_filename.replace("_2", "_1")
        row["outfile2"] = only_or_second_filename
    else:
        row["libtype"] = "unpaired"
        row["fastq_acc"] = only_or_second_filename.rsplit(".fastq.gz")[0]
        row["outfile1"] = only_or_second_filename
    return row

# join acc2ftpDF with metadata df
genome2acc2ftpDF =  m.merge(acc2ftpDF, on = "Sample_accession")

# build fastq_acc and output files
genome2acc2ftpDF= genome2acc2ftpDF.apply(generate_fastq_acc_and_output, axis=1)

# each fastq_acc should have the same associated ftp_links = straight to dict, don't groupby
fastq_acc2ftp = dict(zip(genome2acc2ftpDF.fastq_acc,genome2acc2ftpDF.ftp_links))

# To map genome to multiple input files: 

#first, group by Genome and join all fastq_acc, ignoring libytpe. This will give you _1 files
genome_info1 = genome2acc2ftpDF.groupby("Genome").agg(library1=('fastq_acc','_'.join), r1=('outfile1', list)).reset_index()
# to get _2 files, subset for paired files, then group + aggregate to list.
genome_info2 = genome2acc2ftpDF[genome2acc2ftpDF["libtype"] == "paired"].groupby(["Genome"]).agg(library2=('fastq_acc','_'.join), r2=('outfile2', list)).reset_index()

# now merge these two DF's on the Genome column.
merged_genome_info = genome_info1.merge(genome_info2, on = "Genome", how="outer") # outer takes union, not intersection
# now we have some pesky nan's in the r2, library2 columns from se-only Genomes. fill with ""
merged_genome_info.fillna("", inplace=True)
# if desired you could also print this to a csv, for file provenance metadata storage: 
#merged_genome_info.to_csv(index=False) # and or select cols. you may just want "Genome", "r1", "r2"

# make genome: download file dicts
genome2r1=dict(zip(merged_genome_info.Genome, merged_genome_info.r1))
genome2r2=dict(zip(merged_genome_info.Genome, merged_genome_info.r2))



rule all:
    input:
        #expand("inputs/mgnify_genomes/human-gut/v1.0/{genome}.gff3.gz", genome = GENOMES)
        expand("outputs/charcoal/{genome}.fa.clean.fa", genome = GENOMES),
        "outputs/charcoal_clean_checkm_qa/qa.tsv",
        expand("outputs/charcoal_dirty_prokka/{genome}.fna", genome = GENOMES),
        "outputs/charcoal_dirty_checkm_qa/qa.tsv",
        expand("outputs/magpurify_clean/{genome}_magpurify.fna", genome = GENOMES), 
        expand("outputs/refinem/bams/{genome}.bam", genome = genome2r1.keys())
        

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


rule refinem_download_se_fastq_files:
    output: 
        "inputs/mgnify_raw_reads/se/{fastq_acc}.fastq.gz", 
    wildcard_constraints:
        # no underscore! so avoid issues with _1,_2 
        fasttq_acc='[a-z0-9]+'
    params:
         ftp= lambda w: fastq_acc2ftp[w.fastq_acc][0],
    shell:
        """
        curl -o {output} {params.ftp}
        """

rule refinem_download_fastq_files:
    params:
        ftp_1 = lambda w: fastq_acc2ftp[w.fastq_acc][0],
        ftp_2 = lambda w: fastq_acc2ftp[w.fastq_acc][1]
    output: 
        r1="inputs/mgnify_raw_reads/pe/{fastq_acc}_1.fastq.gz",
        r2="inputs/mgnify_raw_reads/pe/{fastq_acc}_2.fastq.gz",
    wildcard_constraints:
        fastq_acc='[a-z0-9]+'
    shell:
        """
        curl -o {output.r1} {params.ftp_1}
        curl -o {output.r2} {params.ftp_2}
        """ 

def get_cat1_input(w):
    cat_files=[]
    input_files = genome2r1[w.genome]
    for infile in input_files:
        if "_1.fastq.gz" in infile:
            cat_files+=[f"inputs/mgnify_raw_reads/pe/{infile}"]
        else:
            cat_files+=[f"inputs/mgnify_raw_reads/se/{infile}"]
    print(cat_files)
    return cat_files 

rule refinem_cat_libraries_R1:
    input: get_cat1_input
    output: "inputs/cat/{genome}_1.fastq.gz"
    shell:
        """
        cat {input} > {output}
        """

def get_cat2_input(w):
    # look out for "" from nan's!
    cat_files=[]
    input_files = genome2r2[w.genome]
    for infile in input_files:
        if infile:
            cat_files+=[f"inputs/mgnify_raw_reads/pe/{infile}"]
    print(cat_files)
    return cat_files 

rule refinem_cat_libraries_R2:
    input: get_cat2_input
    output: "inputs/cat/{genome}_2.fastq.gz"
    shell:
        """
        cat {input} > {output}
        """

rule refinem_index_genomes:
    input: 
        genome ="outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa",
    output: "outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa.bwt"
    conda: "envs/refinem.yml"
    shell:'''
    bwa index {input}
    '''


def get_align_info(w):
    genome =f"outputs/mgnify_genomes/human-gut/v1.0/{w.genome}.fa"
    indx = f"outputs/mgnify_genomes/human-gut/v1.0/{w.genome}.fa.bwt"
    reads = [f"inputs/cat/{w.genome}_1.fastq.gz"]
    if genome2r2.get(w.genome):
        # it's possible there may only be one read if all are se 
        reads+= [f"inputs/cat/{w.genome}_2.fastq.gz"]
    return {"genome":genome, "indx": indx, "reads": reads }


rule refinem_align_reads:
    input: unpack(get_align_info)
    output: "outputs/refinem/bams/{genome}.bam"
    conda: "envs/bwa.yml"
    benchmark: "benchmarks/refinem_{genome}_bwa_align_reads.txt"
    shell:'''
    bwa mem {input.genome} {{" ".join(reads)}} | samtools sort -o {output} -")
    ''' 

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
    refinem outliers {input} {params.outdir}
    '''

rule refinem_filter_bins_stats:
    input: "outputs/refinem/outliers/outliers.tsv"
    output: directory("outputs/refinem/filter_bin_stats/") # need to decide if I want all bins here or just modified bins. if mod, probs directory().
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        outdir = "outputs/refinem/filter_bins_stats/"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_filter_bins_stats.txt"
    shell:'''
    refinem filter_bins {params.bindir} {input} {params.outdir} --modified_only # flag makes this rule only output modified bins
    '''

rule refinem_call_genes:
    input: 
        expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
    output: directory("outputs/refinem/call_genes") # NEED TO FILL IN OUTPUT FILE (run first, to determine name of an output file)
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        genedir = "outputs/refinem/call_genes/"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_call_genes.txt"
    shell:'''
    refinem call_genes -c 4 {params.bindir} {params.genedir}
    '''

rule refinem_download_protein_db:
    output:"inputs/refinem_db/gtdb_r89_protein_db.2019-09-27.faa.gz"
    shell:'''
    wget -O {output} https://data.ace.uq.edu.au/public/misc_downloads/refinem/gtdb_r89_protein_db.2019-09-27.faa.gz
    '''

rule refinem_download_taxonmy_db:
    output: "inputs/refinem_db/gtdb_r89_taxonomy.2019-09-27.tsv"
    shell:'''
    wget -O {output} https://data.ace.uq.edu.au/public/misc_downloads/refinem/gtdb_r89_taxonomy.2019-09-27.tsv
    '''

rule refinem_taxon_profile:
    input: 
        genedir = directory("outputs/refinem/call_genes"),
        scaffold_stats = "outputs/refinem/scaffold_stats/scaffold_stats.tsv",
        db = "inputs/refinem_db/gtdb_r89_protein_db.2019-09-27.faa.gz",
        tax = "inputs/refinem_db/gtdb_r89_taxonomy.2019-09-27.tsv"
    output: directory("outputs/refinem/taxon_profile/") # NEED TO ADD A SPECIFIC FILE NAME
    params:
        genedir = "outputs/refinem/call_genes",
        outdir = "outputs/refinem/taxon_profile"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_taxon_profile.txt"
    shell:'''
    refinem taxon_profile -c 4 {params.genedir} {input.scaffold_stats} {input.db} {input.tax} {params.outdir}
    '''

rule refinem_taxon_filter:
    input: directory("outputs/refinem/taxon_profile/") # NEED TO ADD A SPECIFIC FILE NAME TO MATCH OUTPUT OF ABOVE RULE 
    output: "outputs/refinem/taxon_profile/taxon_filter.tsv" 
    params: indir = "outputs/refinem/taxon_profile"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_taxon_filter.txt"
    shell:'''
    refinem taxon_filter -c 4 {params.indir} {output}
    '''

rule refinem_filter_bins_taxon:
    input: 
        tax_filt = "outputs/refinem/taxon_profile/taxon_filter.tsv",
        genome = expand("outputs/mgnify_genomes/human-gut/v1.0/{genome}.fa", genome = GENOMES),
    output: directory("outputs/refinem/filter_bins_taxon") # NEED TO ADD SPECIFIC FILE OUTPUTS
    params: 
        bindir = "outputs/mgnify_genomes/human-gut/v1.0",
        outdir = "outputs/refinem/filter_bins/taxon"
    conda: "envs/refinem.yml"
    benchmark: "benchmarks/refinem_filter_bins_taxon.txt"
    shell: '''
    refinem filter_bins {params.bindir} {input.tax_filt} {params.outdir}
    '''
