# 2020 Charcoal Paper

Analyses and writeup of [charcoal](github.com/dib-lab/charcoal).

## Getting started

To reexecute the analyses in this repository:

```
conda env create -f environment.yml
conda activate charcoal_paper

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=300000 --cluster "sbatch -t 20160 -J char -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k
```

Note that due to long DAG solve times for running charcoal on all genomes in GTDB rs202, we split the charcoal runs into twenty subsets. 
Configuration files for each run can be found in `inputs/charcoal_conf`, while run instructions can be found in `run_charcoal_wo_smk.sh`.

## Outline 

Rough results outline, computation for which is accomplished in this repository.

1. GTDB rs202 decontamination at order rank
   - charcoal assigns accurate lineage 
   - number of genomes contaminated
   - average % of genome contaminated, average contaminant bps
   - distribution of contamination 
       - MAGs vs isolates vs RefSeq 
       - GTDB rep vs. non rep
       - by taxonomy
       - pairs of contaminants occur more than expected by chance
   - comparison to checkm contamination
   - number of shared 31-mers decreases after GTDB decontamination 
2. Run charcoal on sets of MAGs from different environments
   - mgnify
   - TARA oceans
3. Use case of Eukaryotic decontamination
