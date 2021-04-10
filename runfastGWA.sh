#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=GCTA-fastGWA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=167:00:00
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
grmdir=("/scratch/mf91122/CCC/exomeQC200k/GCTA-GRM-CALL/sparseGRM")
outdir=("/scratch/mf91122/CCC/exomeQC200k/G.fastGWA-04082021")

phenotypes=("Creatinine_resinv" "Cystatin_C_resinv" "eGFR1_resinv" "eGFR2_resinv" "eGFR3_resinv" \ 
              "eGFR4_resinv" "BMI_resinv" "SBP_resinv" "DBP_resinv" "HbA1c_resinv" "LDL_resinv" "HDL_resinv" \
              "TAGs_resinv" "TC_resinv" "Waist_circumference_resinv")


#phenotypes=("Creatinine_resinv")

for j in ${phenotypes[@]} 
	do

mkdir -p $outdir/$j

$gctadir/gcta64 \
--mbfile $genoindir/geno_chrs.txt \
--grm-sparse $grmdir/sparse_grm \
--fastGWA-mlm \
--pheno $phenodir/$j.phen \
--qcovar $phenodir/KidneyPhenoCovarPC.txt \
--threads 8 \
--out $outdir/$j/geno_assoc1

done
