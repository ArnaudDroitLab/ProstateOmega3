for i in output/*/STARAligned.out.bam
do
    sample=`echo $i | sed -e 's/output\/\(.*\)\/STARAligned.out.bam/\1/'`
    script=output/$sample/filter_count.sh
    cat <<'EOF' > $script
#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

module load samtools
module load python

samtools view -b -F 4 -q 20 output/$sample/STARAligned.out.bam > output/$sample/STARAligned.filtered.bam
samtools sort -@ 8 output/$sample/STARAligned.filtered.bam > output/$sample/STARAligned.filtered.sorted.bam
samtools index -@ 8 output/$sample/STARAligned.filtered.bam
#rm output/$sample/STARAligned.filtered.bam
htseq-count --order pos -s reverse -f bam output/$sample/STARAligned.filtered.bam output/introns_exons.gtf > output/$sample/htseq_count_intron.txt
EOF

    sbatch --export=sample=$sample --job-name=$sample.htseq -o $script.stdout -e $script.stderr -D `pwd` $script 
done

