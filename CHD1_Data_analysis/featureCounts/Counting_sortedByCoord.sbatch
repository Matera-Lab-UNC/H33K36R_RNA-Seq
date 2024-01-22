#!/bin/bash

#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts_%A.out
#SBATCH --error=featureCounts_%A.err
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

cd /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/Bam

ml subread/2.0.3

bam1=SRR11237854_Aligned.sortedByCoord.out.bam
bam2=SRR11237855_Aligned.sortedByCoord.out.bam
bam3=SRR11237856_Aligned.sortedByCoord.out.bam
bam4=SRR11237857_Aligned.sortedByCoord.out.bam
bam5=SRR11237858_Aligned.sortedByCoord.out.bam
bam6=SRR11237859_Aligned.sortedByCoord.out.bam
bam7=SRR11237860_Aligned.sortedByCoord.out.bam
bam8=SRR11237861_Aligned.sortedByCoord.out.bam
bam9=SRR11237862_Aligned.sortedByCoord.out.bam

gtf=/proj/materalb/users/John/STAR_Genome_Index_Dm6.38/dm6_fb.filtered.gtf

featureCounts \
-a $gtf \
-o /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/CHD1_Data_analysis/featureCounts/sortedByCoord_counts.txt \
-T 10 \
-g gene_id \
-s 0 \
$bam1 $bam2 $bam3 $bam4 $bam5 $bam6 $bam7 $bam8 $bam9
