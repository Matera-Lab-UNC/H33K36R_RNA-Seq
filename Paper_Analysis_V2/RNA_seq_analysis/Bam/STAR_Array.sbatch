#!/bin/bash

#SBATCH --job-name=STAR_Align
#SBATCH --output=STAR_Align_%A_%a.out
#SBATCH --ntasks=16
#SBATCH --mem=40G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

cd /proj/materalb/users/John/Merged_lane_files

module load star/2.7.7a

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./fastq.list) 

STAR --quantMode TranscriptomeSAM GeneCounts \
--genomeDir /proj/materalb/bin/Indexes/STAR_Indexes/GenomeDir_dm6.38_indices_for_alignment \
--readFilesIn $(ls | grep -e ${line} | grep -e "R1.fastq.gz") $(ls | grep -e ${line} | grep -e "R2.fastq.gz") \
--outFileNamePrefix /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/Paper_Analysis_V2/RNA_seq_analysis/Bam/${line}_ \
--readFilesCommand zcat \
--sjdbGTFfile /proj/materalb//bin/Indexes/STAR_Indexes/GenomeDir_dm6.38_indices_for_alignment/dm6_fb.filtered.gtf \
--runThreadN 16 \
--outFilterMismatchNmax 2 \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate
