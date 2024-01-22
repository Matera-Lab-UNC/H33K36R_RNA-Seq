#!/bin/bash

#SBATCH --job-name=VDBvalidate
#SBATCH --output=VDBvalidate_%A.out
#SBATCH --time=05:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bdmcmi@ad.unc.edu

module load sratoolkit/3.0.1

vdb-validate $(<SRR_Acc_List.txt)
