#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=BOLTtest2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=167:00:00
#SBATCH --mem=1000000
#SBATCH --output=BOLT.%j.out
#SBATCH --error=BOLT.%j.err
#SBATCH --constraint=Intel

cd /work/kylab/mike/CCC/BOLT-03012020
#https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4_manual.pdf

genoindir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES/bgenKEEP/bfile/combine/bgen")
phenodir=("/scratch/mf91122/CCC/pheno200kresid")
outdir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES/BOLT3-1mil")
boltdir=("/home/mf91122/BOLT-LMM/BOLT-LMM_v2.3.4")

#j="Creatinine_resinv"
phenotypes=("Cystatin_C_resinv" "eGFR1_resinv" "eGFR2_resinv" "eGFR3_resinv" "eGFR4_resinv")

for j in ${phenotypes[@]} 
	do

mkdir -p $outdir/$j

$boltdir/bolt \
--bfile=/scratch/mf91122/CCC/exomeQC200k/Prune1/bfileKEEP/combine/merged \
--phenoFile="$phenodir"/KidneyPhenoFullpc.txt \
--phenoCol="$j" \
--covarFile="$phenodir"/KidneyPhenoFullpc.txt \
--qCovarCol=PC{1:10} \
--bgenFile="$genoindir"/combine.bgen \
--sampleFile="$genoindir"/combine.sample \
--lmm \
--numThreads=16 \
--LDscoresFile="$boltdir"/tables/LDSCORE.1000G_EUR.tab.gz \
--bgenMinINFO=0.3 \
--statsFile="$outdir"/"$j"/BOLT1-statsFile \
--statsFileBgenSnps="$outdir"/"$j"/BOLT1-statsFile-BgenSnps

done
