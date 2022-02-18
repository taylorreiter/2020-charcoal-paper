srun -p bmh -J char600 -t 600:00:00 --mem=16gb -c 2 --pty bash
conda activate charcoal_paper
conda activate /home/tereiter/github/2022-microberna/.snakemake/conda/1e7ae02e473e49aa270ffd4484ffd236
export TMPDIR='/scratch/tereiter'
python -c "import os; print(os.environ['TMPDIR'])"
python -m charcoal run inputs/charcoal_conf/gtdb_rs202_genomes1.conf -j 200 clean --nolock --use-conda --latency-wait 15 --rerun-incomplete -k --cluster "sbatch -t 90:00 -J char -p bml -n 1 -N 1 -c 1 --mem=34Gb"
