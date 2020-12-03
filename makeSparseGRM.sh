#!/bin/bash
#SBATCH --partition=highmem_30d_p
#SBATCH --job-name=GCTA-sparseGRM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=600:00:00
#SBATCH --mem=600000
#SBATCH --output=GCTA-sparse1.%j.out
#SBATCH --error=GCTA-sparse1.%j.err

#use a new GRM function in GCTA using the BLAS packages available in the Intel MKL library to compute a GRM
#https://cnsgenomics.com/software/gcta/#MakingaGRM
#https://cnsgenomics.com/software/gcta/#fastGWA

cd /work/kylab/mike/CCC/C2.GCTA-GRM

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
genoindir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES")

outdir=("/scratch/mf91122/CCC/exomeQC200k/C2.GCTA-GRM/sparse")

mkdir -p $outdir

$gctadir/gcta64 \
--grm $outdir/../geno_grm \
--make-bK-sparse 0.05 \
--out $outdir/sp_grm
