for i in output/*/STARAligned.out.bam
do
    sample=`echo $i | sed -e 's/output\/\(.*\)\/STARAligned.out.bam/\1/'`
    script=output/$sample/count.sh
    cat <<'EOF' > $script
#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

module load python

htseq-count --order pos -s reverse -f bam output/$sample/STARAligned.filtered.bam output/introns_combined.gtf > output/$sample/htseq_count_intron.txt
EOF

    sbatch --export=sample=$sample --job-name=$sample.htseq -o $script.stdout -e $script.stderr -D `pwd` $script 
done

