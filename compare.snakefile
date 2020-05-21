

rule all:
    input:
        


rule run_charcoal:
    input:
        mag_list = ,
        conf = ,
    output: "{}"
    conda: "envs/charcoal.yml"
    benchmark: "benchmarks/charcoal.benchmark"
    shell:'''
    charcoal run {input.conf} -j 8
    '''

rule run_checkm_lineage_qf_clean:
    """
    Reports general summary statistics (% complete, % contamination) for genome (bins).
    """
    input: "outputs/charcoal/{}.clean.fa.gz"
    output: "outputs/charcoal_clean_checkm/completeness.tsv"
    params: 
        indir = "outputs/charcoal",
        threads = "4"
    conda: "envs/checkm.yml"
    benchmark: "benchmarks/{}_checkm_clean.benchmark"
    shell:'''
    checkm lineage_wf \
        --file {output} \
        --tab_table \
        --extension clean.fa.gz \
        --threads {params.threads} \
        {params.indir} 
    '''     

rule run_checkm_qa_clean:
    """
    Reports name of contig on which contaminant marker gene resides, as well as marker gene identity.
    """
    input: marker_file = "charcoal_clean_checkm/lineage.ms"
    output:
    conda: "envs/checkm.yml"
    benchmark: "benchmarks/{}_checkm_qa_clean.benchmark"
    params: 
        threads = "4",
        indir = "charcoal_clean_checkm" 
    shell:'''
    checkm qa -o 8 -f {output} --tab_table -t {params.threads} {input.marker_file} {params.indir}
    '''

rule run_prokka_dirty:
    input: "outputs/charcoal/{}.dirty.fa.gz"
    output: "outputs/charcoal_dirty_prokka/{}"
    params: outdir = "outputs/charcoal_dirty_prokka" 
    conda: "envs/prokka.yml"
    benchmark: "benchmarks/{}_prokka_dirty.benchmark"
    shell:'''
    prokka {input} --outdir {params.outdir} --prefix {wildcards.*} --metagenome --force --locustag {wildcards.*}
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
