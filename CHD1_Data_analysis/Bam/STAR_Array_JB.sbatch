#!/bin/bash

#SBATCH --job-name=STAR_align
#SBATCH --output=STAR_align_%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=16
#SBATCH --array=1-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

cd /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Fastq/Trim_Fastq

module load star/2.7.7a

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./SRR_Acc_List.txt)

STAR --quantMode TranscriptomeSAM GeneCounts \
--genomeDir /proj/materalb/bin/Indexes/STAR_Indexes/GenomeDir_dm6.38_indices_for_alignment \
--readFilesIn $(ls | grep -e ${line} | grep -e "_trim.fastq.gz") \
--outFileNamePrefix /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Bam/${line}_ \
--readFilesCommand zcat \
--sjdbGTFfile /proj/materalb/users/John/STAR_Genome_Index_Dm6.38/dm6_fb.filtered.gtf \
--runThreadN 16 \
--outFilterMismatchNmax 2 \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate
