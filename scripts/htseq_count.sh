#!/bin/bash

module load python

htseq-count --order pos -s reverse -f bam output/P020_V06_BP2_S75/STARAligned.filtered.bam ref/Homo_sapiens.GRCh38.94.chr.gtf > output/P020_V06_BP2_S75/htseq_count_ensembl-reverse.txt
htseq-count --order pos -s reverse -f bam output/P020_V06_BP2_S75/STARAligned.filtered.bam output/P020_V06_BP2_S75/stringtie.gtf > output/P020_V06_BP2_S75/htseq_count_stringtie-reverse.txt

htseq-count --order pos -s reverse -f bam output/P020_V06_BP2_S75/STARAligned.filtered.bam combined.gtf > output/P020_V06_BP2_S75/htseq_count_intron.txt
