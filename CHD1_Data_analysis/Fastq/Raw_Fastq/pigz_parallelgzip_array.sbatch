#!/bin/bash

#SBATCH --job-name=pigz
#SBATCH --output=pigz_%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=10
#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

module load pigz/2.7

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./SRR_Acc_List.txt)

pigz -k ${line}.fastq
