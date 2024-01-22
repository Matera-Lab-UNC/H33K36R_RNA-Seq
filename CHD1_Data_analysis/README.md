# Directory created on 06-26-2023
## By Benjamin McMichael

------------------------------------------------------------------------------------------------------------------
# Contents:

## Fastq files loaded from SRA using SRA-Toolkit.
Files are from:
'''
	Schoberleitner I, Bauer I, Huang A, Andreyeva EN, Sebald J, Pascher K, Rieder D, Brunner M, Podhraski V, Oemer G, Cázarez-García D, Rieder L, Keller MA, Winkler R, Fyodorov DV, Lusser A. CHD1 controls H3.3 incorporation in adult brain chromatin to maintain metabolic homeostasis and normal lifespan. Cell Rep. 2021 Oct 5;37(1):109769. doi: 10.1016/j.celrep.2021.109769. PMID: 34610319; PMCID: PMC8607513.
 '''

The prefetch command was used to load SRA data for all samples.
The fasterq-dump command was used to load fastq files for each sample.
		Prefetch should be performed before fasterq-dump to speed up file downloads.
Fastq files were zipped using pigz.
		Parallel implementation of gzip.
		fasterq-dump does NOT have --gzip option.

Fastq files were quality AND adapter trimmed using bbduk.

Bam files were aligned using STAR.
	Used the same STAR parameters John used for his fly head files.
	Confirmed this by checking his analysis in matera lab project directory.
	Aligned to the dm6 index, and used the filtered dm6 GTF file generated by Harmony. 