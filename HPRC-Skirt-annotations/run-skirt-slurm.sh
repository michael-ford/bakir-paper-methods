#! /bin/bash
#SBATCH --job-name=run-skirt
#SBATCH --cpus-per-task=5
#SBATCH --mem=4G
#SBATCH --gres=lscratch:200
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL

module load nextflow
export TMPDIR=/lscratch/$SLURM_JOB_ID

nextflow run run-skirt.nf \
-profile biowulf \
