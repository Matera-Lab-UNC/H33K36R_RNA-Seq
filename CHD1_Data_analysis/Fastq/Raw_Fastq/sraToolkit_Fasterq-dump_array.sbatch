#!/bin/bash

#SBATCH --job-name=fasterq-dump
#SBATCH --output=fasterq-dump_%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=6
#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

module load sratoolkit/3.0.1

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./SRR_Acc_List.txt)

fasterq-dump --progress --threads 6 ${line}
