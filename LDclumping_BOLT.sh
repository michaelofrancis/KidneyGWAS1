#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=LD_clumping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=LD_clumping.%j.out
#SBATCH --error=LD_clumping.%j.err

cd /work/kylab/mike/CCC/LD_clumping

ml PLINK/1.9b_5-x86_64


indir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES/BOLT3")
outdir=("/scratch/mf91122/CCC/exomeQC200k/BOLT-LD-clumping")

mkdir -p $outdir

#Clumping analysis criteria: P-value threshold = 0.05 divided by the number of variants tested for each trait 
#Number of variants tested for each trait:
##wc -l BOLT1-statsFile-BgenSnps ##2865352
#0.05/ 2865352 = 1.74e-08
#window size = 5 Mb
#LD r2 threshold = 0.01. 


plink \
--clump $indir/testLDinput \
--clump-p1 1.74e-08 \
--clump-r2 0.01 \
--clump-kb 5000 \
--out "$outdir"/Creatinine-LDclump
