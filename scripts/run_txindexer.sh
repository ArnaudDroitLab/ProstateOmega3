#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

export PATH=/home/efournie/FuSeq_v1.1.0_linux_x86-64/linux/bin:$PATH

TxIndexer -p 8 -t ref/Homo_sapiens.GRCh38.cdna.all.nover.fa -o ref/FuSeq_idx
