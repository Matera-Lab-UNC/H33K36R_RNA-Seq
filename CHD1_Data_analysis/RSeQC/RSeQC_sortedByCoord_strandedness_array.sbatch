#!/bin/bash

#SBATCH --job-name=RSeQC_strandedness
#SBATCH --output=RSeQC_strandedness_%A_%a.out
#SBATCH --error=RSeQC_strandedness_%A_%a.err
#SBATCH --time=30:00
#SBATCH --mem=10G
#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

ml rseqc/4.0.0

cd /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Bam

declare -a samples=(ls *.sortedByCoord.out.bam)

c=${samples[$SLURM_ARRAY_TASK_ID]}

infer_experiment.py \
--input-file=${c} \
--refgene=/nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/RSeQC/dm6_fb.filtered.bed
