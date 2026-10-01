#! /bin/bash -login
#SBATCH -J phaserPopCisVar
#SBATCH -t 4320
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8gb
#SBATCH -o "%x-%j.out"

## Exit if any command fails
set -e

## Create and activate virtual environment with requirements
conda activate /gpfs/commons/home/mbond/.conda/snakemake7_env

## Execute asePipe snakemake workflow
slurm_config_folder="config/snakeConfig"
snakemake --snakefile workflow/cisVar.smk --profile $slurm_config_folder --nolock -k --rerun-triggers mtime --rerun-incomplete

## Success message
echo "Entire workflow completed successfully!"
