#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=GCTA-sparseGRM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=72:00:00
#SBATCH --mem=600000
#SBATCH --output=GCTA-sparse1.%j.out
#SBATCH --error=GCTA-sparse1.%j.err

#use a new GRM function in GCTA using the BLAS packages available in the Intel MKL library to compute a GRM
#https://cnsgenomics.com/software/gcta/#MakingaGRM
#https://cnsgenomics.com/software/gcta/#fastGWA

cd /work/kylab/mike/CCC/GCTA-GRM-CALL

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")
genoindir=("/scratch/mf91122/CCC/UKB_call")

outdir=("/scratch/mf91122/CCC/exomeQC200k/GCTA-GRM-CALL")

mkdir -p $outdir

$gctadir/gcta64 \
--mbfile /scratch/mf91122/CCC/UKB_call/chrs.txt \
--make-grm-part 3 1 \
--extract /scratch/mf91122/CCC/exomeQC200k/Prune3/Prune3.in.txt \
--thread-num 16 \
--out $outdir/grm

$gctadir/gcta64 \
--mbfile /scratch/mf91122/CCC/UKB_call/chrs.txt \
--make-grm-part 3 2 \
--extract /scratch/mf91122/CCC/exomeQC200k/Prune3/Prune3.in.txt \
--thread-num 16 \
--out $outdir/grm

$gctadir/gcta64 \
--mbfile /scratch/mf91122/CCC/UKB_call/chrs.txt \
--make-grm-part 3 3 \
--extract /scratch/mf91122/CCC/exomeQC200k/Prune3/Prune3.in.txt \
--thread-num 16 \
--out $outdir/grm

# Merge all the parts together (Linux, Mac)
#cat test.part_3_*.grm.id > test.grm.id
#cat test.part_3_*.grm.bin > test.grm.bin
#cat test.part_3_*.grm.N.bin > test.grm.N.bin


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-Second script =-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=GCTA-sparseGRM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=72:00:00
#SBATCH --mem=100000
#SBATCH --output=GCTA-sparse1.%j.out
#SBATCH --error=GCTA-sparse1.%j.err

#use a new GRM function in GCTA using the BLAS packages available in the Intel MKL library to compute a GRM
#https://cnsgenomics.com/software/gcta/#MakingaGRM
#https://cnsgenomics.com/software/gcta/#fastGWA

cd /work/kylab/mike/CCC/GCTA-GRM-CALL

gctadir=("/home/mf91122/GCTA/gcta_1.93.2beta")

outdir=("/scratch/mf91122/CCC/exomeQC200k/GCTA-GRM-CALL/sparseGRM")

mkdir -p $outdir

$gctadir/gcta64 \
--grm $outdir/../fullGRM/grm \
--make-bK-sparse 0.05 \
--out $outdir/sparse_grm
