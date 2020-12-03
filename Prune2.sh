#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=Prune2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=Prune2.%j.out
#SBATCH --error=Prune2.%j.err
#SBATCH --array=11-22

cd /work/kylab/mike/CCC/C.QC
#ml PLINK/2.00-alpha2.3-x86_64-20200914-dev
ml PLINK/1.9b_5-x86_64

genoindir=("/scratch/mf91122/CCC/UKB_call")
outdir=("/scratch/mf91122/CCC/exomeQC200k/Prune2")

mkdir -p $outdir

i=$SLURM_ARRAY_TASK_ID

plink \
--bed "$genoindir"/ukb_cal_chr"$i"_v2.bed \
--fam "$genoindir"/famfiles/ukb48818_cal_chr"$i"_v2_s488282.fam \
--bim "$genoindir"/bimfiles/ukb_snp_chr"$i"_v2.bim \
--geno 0.01 \
--hwe 1e-06 \
--maf 0.4 \
--exclude range Prune1excludeset.txt \
--keep /scratch/mf91122/CCC/pheno200kresid/KidneyPhenoIDFullFinal.txt \
--indep-pairwise 1000 50 0.05 \
--out "$outdir"/chr"$i"
