#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=QC_UKB_WES
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=QC_UKB_WES.%j.out
#SBATCH --error=QC_UKB_WES.%j.err

cd /work/kylab/mike/CCC/C.QC
#ml PLINK/2.00-alpha2.3-x86_64-20200914-dev
ml PLINK/1.9b_5-x86_64

genoindir=("/scratch/mf91122/UKBimputation/exome_sequences/200k")
outdir=("/scratch/mf91122/CCC/exomeQC200k/QC_UKB_WES/")

chr=(22 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
#chr=(22)


for i in ${chr[*]}
        do

plink \
--bed "$genoindir"/ukb23155_c"$i"_b0_v1.bed \
--fam "$genoindir"/ukb23155_c1_b0_v1_s200631.fam \
--bim "$genoindir"/bim/UKBexomeOQFE_chr"$i".bim \
--mind 0.05 \
--geno 0.05 \
--hwe 1e-06 \
--mac 10 \
--make-bed \
--out "$outdir"/chr"$i"

done #end chromosome loop
