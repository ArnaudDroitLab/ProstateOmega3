#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

mkdir -p ref/STAR_index

module load star/2.6.1a

STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ref/STAR_index \
    --genomeFastaFiles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa input/ERCC.fa ref/HS_rRNA.fa \
    --sjdbGTFfile ref/Homo_sapiens.GRCh38.94.chr.gtf \
    --sjdbOverhang 124