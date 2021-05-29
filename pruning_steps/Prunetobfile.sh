#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=Prunetobfile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=Prunetobfile.%j.out
#SBATCH --error=Prunetobfile.%j.err
#SBATCH --array=1-22


cd /work/kylab/mike/T-WES-GWAS/Prune
#ml PLINK/2.00-alpha2.3-x86_64-20200914-dev
ml PLINK/1.9b_5-x86_64

genoindir=("/scratch/mf91122/CCC/UKB_call")
outdir=("/scratch/mf91122/T-WES-GWAS")

i=$SLURM_ARRAY_TASK_ID

for j in {1..3}
	do

mkdir -p $outdir/Prune"$j"/bfileKEEP/combine

#Make by chromosome
plink \
--bfile /scratch/mf91122/CCC/UKB_call/chr"$i" \
--extract /scratch/mf91122/T-WES-GWAS/Prune"$j"/chr"$i".prune.in \
--keep /scratch/mf91122/T-WES-GWAS/phenotables/T-WES-IDmodel0_05292021.txt \
--make-bed \
--out "$outdir"/Prune"$j"/bfileKEEP/chr"$i"


#Merge
plink \
--make-bed \
--merge-list "$outdir"/Prune"$j"/bfileKEEP/list.txt \
--out "$outdir"/Prune"$j"/bfileKEEP/combine/merged

done #end j loop
