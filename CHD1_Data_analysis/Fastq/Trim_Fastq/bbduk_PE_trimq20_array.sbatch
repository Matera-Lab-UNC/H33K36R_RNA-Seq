#!/bin/bash

#SBATCH --job-name=bbduk
#SBATCH --output=bbduk_%A_%a.out
#SBATCH --error=bbduk_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=25G
#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

cd /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Fastq/Raw_Fastq

module load bbmap/39.01

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./SRR_Acc_List.txt)

bbduk.sh \
in=${line}.fastq.gz \
out1=/nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Fastq/Trim_Fastq/${line}_trim.fastq.gz \
qtrim=rl \
trimq=20 \
ref=adapters \
tbo \
tpe
