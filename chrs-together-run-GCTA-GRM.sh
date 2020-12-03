#!/bin/bash
#SBATCH --partition=highmem_30d_p
#SBATCH --job-name=GCTA-GRM-2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=600:00:00
#SBATCH --mem=600000
#SBATCH --output=GCTA-GRM2.%j.out
#SBATCH --error=GCTA-GRM2.%j.err
#SBATCH --array=18

#use a new GRM function in GCTA using the BLAS packages available in the Intel MKL library to compute a GRM
#https://cnsgenomics.com/software/gcta/#MakingaGRM

cd /work/kylab/mike/CCC/C2.GCTA-GRM


i=$SLURM_ARRAY_TASK_ID

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
genoindir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES")

outdir=("/scratch/mf91122/CCC/exomeQC200k/C2.GCTA-GRM")

mkdir -p $outdir

$gctadir/gcta64 \
--mbfile $genoindir/geno_chrs.txt \
--make-grm \
--thread-num 16 \
--out $outdir/geno_grm
