#!/bin/bash

#SBATCH --job-name=featureCounts_forStrand
#SBATCH --output=featureCounts_forStrand_%A.out
#SBATCH --error=featureCounts_forStrand_%A.err
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

cd /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/Paper_Analysis_V2/RNA_seq_analysis/Bam

ml subread/2.0.3

bam1=MATA_LB_JB_70-y_Adel_f_ACGTCGTT_S1_Aligned.sortedByCoord.out.bam
bam2=MATA_LB_JB_71-y_Adel_m_GTCCTGTT_S2_Aligned.sortedByCoord.out.bam
bam3=MATA_LB_JB_72-o_Adel_f_AGAAGCCT_S3_Aligned.sortedByCoord.out.bam
bam4=MATA_LB_JB_73-o_Adel_m_GAAGATCC_S4_Aligned.sortedByCoord.out.bam
bam5=MATA_LB_JB_74-y_K36_f_TAGCTGAG_S5_Aligned.sortedByCoord.out.bam
bam6=MATA_LB_JB_75-y_K36_m_ACGTCCAA_S6_Aligned.sortedByCoord.out.bam
bam7=MATA_LB_JB_76-o_K36_f_CACACATC_S7_Aligned.sortedByCoord.out.bam
bam8=MATA_LB_JB_77-o_K36_m_CGGATCAA_S8_Aligned.sortedByCoord.out.bam
bam9=MATA_LB_JB_78-y_Adel_f_TCAGCCTT_S9_Aligned.sortedByCoord.out.bam
bam10=MATA_LB_JB_79-y_Adel_m_AAGGCTCT_S10_Aligned.sortedByCoord.out.bam
bam11=MATA_LB_JB_80-o_Adel_f_TGTTCCGT_S11_Aligned.sortedByCoord.out.bam
bam12=MATA_LB_JB_81-o_Adel_m_GGAATGTC_S12_Aligned.sortedByCoord.out.bam
bam13=MATA_LB_JB_82-y_K36_f_CATCCAAG_S13_Aligned.sortedByCoord.out.bam
bam14=MATA_LB_JB_83-y_K36_m_GTCAACAG_S14_Aligned.sortedByCoord.out.bam
bam15=MATA_LB_JB_84-o_K36_f_TCGCTATC_S15_Aligned.sortedByCoord.out.bam
bam16=MATA_LB_JB_85-o_K36_m_AGCCTATC_S16_Aligned.sortedByCoord.out.bam
bam17=MATA_LB_JB_86-y_Adel_f_TCGGATTC_S17_Aligned.sortedByCoord.out.bam
bam18=MATA_LB_JB_87-y_Adel_m_CGGAGTAT_S18_Aligned.sortedByCoord.out.bam
bam19=MATA_LB_JB_88-o_Adel_f_GAACCTTC_S19_Aligned.sortedByCoord.out.bam
bam20=MATA_LB_JB_89-o_Adel_m_AGAGGATG_S20_Aligned.sortedByCoord.out.bam
bam21=MATA_LB_JB_90-y_K36_f_ACGCTTCT_S21_Aligned.sortedByCoord.out.bam
bam22=MATA_LB_JB_91-y_K36_m_CACAGGAA_S22_Aligned.sortedByCoord.out.bam
bam23=MATA_LB_JB_92-o_K36_f_ACGAATCC_S23_Aligned.sortedByCoord.out.bam
bam24=MATA_LB_JB_93-o_K36_m_CCTTCCAT_S24_Aligned.sortedByCoord.out.bam

gtf=/proj/materalb//bin/Indexes/STAR_Indexes/GenomeDir_dm6.38_indices_for_alignment/dm6_fb.filtered.gtf

featureCounts \
-a $gtf \
-o /nas/longleaf/home/bdmcmi/work/Benjamin/flyHead_analysis/Paper_Analysis_V2/RNA_seq_analysis/Counts/featureCounts_sortedByCoord_forStrand.txt \
-T 10 \
-g gene_id \
-p \
--countReadPairs \
-s 1 \
$bam1 $bam2 $bam3 $bam4 $bam5 $bam6 $bam7 $bam8 $bam9 $bam10 $bam11 $bam12 \
$bam13 $bam14 $bam15 $bam16 $bam17 $bam18 $bam19 $bam20 $bam21 $bam22 $bam23 $bam24
