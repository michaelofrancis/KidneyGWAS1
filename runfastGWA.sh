#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=GCTA-fastGWA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=500000
#SBATCH --output=GCTA-fastGWA.%j.out
#SBATCH --error=GCTA-fastGWA.%j.err

#fastGWA (https://cnsgenomics.com/software/gcta/#fastGWA) is a new method that can 
#quickly finish GWAS using mixed linear regression model for large cohort, while keeping 
#low false positive rate, controlling for population stratification by PC and relatedness 
#by a sparse genetic relationship matrix.


cd /work/kylab/mike/CCC/G.fastGWA

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
genoindir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES")
phenodir=("/scratch/mf91122/CCC/pheno200kresid/phen")
grmdir=("/scratch/mf91122/CCC/exomeQC200k/C2.GCTA-GRM/sparse")
outdir=("/scratch/mf91122/CCC/exomeQC200k/G.fastGWA")


mkdir -p $outdir

$gctadir/gcta64 \
--mbfile $genoindir/geno_chrs.txt \
--grm-sparse $grmdir/sp_grm \
--fastGWA-mlm \
--pheno $phenodir/Creatinine_resinv.phen \
--qcovar $phenodir/KidneyPhenoCovarPC.txt \
--threads 16 \
--out $outdir/geno_assoc1
