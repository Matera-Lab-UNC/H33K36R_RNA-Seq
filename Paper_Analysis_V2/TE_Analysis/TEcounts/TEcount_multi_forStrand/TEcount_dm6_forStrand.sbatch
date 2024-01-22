#!/bin/bash

#SBATCH --job-name=TE_count_forStrand
#SBATCH --output=TE_count_forStrand_%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --array=1-24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

#cd /work/users/b/d/bdmcmi/Benjamin/flyHead_analysis/Paper_Analysis_V2/RNA_seq_analysis/Bam

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./bam.list)

singularity exec --bind /work/users/b/d/bdmcmi/Benjamin/flyHead_analysis/Paper_Analysis_V2/TE_Analysis/TEcounts/TEcount_multi_forStrand tetranscripts.sif TEcount \
  --format BAM \
  -b ${line}.sortedByCoord.out.bam \
  --GTF /work/users/b/d/bdmcmi/Benjamin/flyHead_analysis/Paper_Analysis_V2/TE_Analysis/TEcounts/TEcount_multi_forStrand/dm6_fb.filtered.gtf \
  --TE /work/users/b/d/bdmcmi/Benjamin/flyHead_analysis/Paper_Analysis_V2/TE_Analysis/TEcounts/TEcount_multi_forStrand/dm6_rmsk_TE_rmCHRprefix.gtf \
  --stranded forward \
  --sortByPos \
  --project TEcount_${line} \
\
