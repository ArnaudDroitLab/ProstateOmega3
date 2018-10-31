module load trimmomatic/0.36
module load kallisto/0.44.0

# Make a list of all unique samples.
rm -f tmp_file_list
curdir=`pwd`
for i in raw/180925_D00487S_0161_ACCR9DANXX/ raw/180925_D00487S_0162_BCCU7UANXX/ raw/181003_D00487S_0163_ACCU7FANXX/
do
    pushd $i
    ls * >> $curdir/tmp_file_list
    popd
done

# Make a list of all samples.
sort < tmp_file_list | perl -ne 'print "$1\n" if /(.+)_L..._R._001.fastq.gz/' | sort | uniq > unique_sample_list

while read sample; do
    mkdir -p output/$sample
    script=output/$sample/trim_kallisto.sh
    cat <<'EOF' > $script
#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

module load trimmomatic/0.36
module load kallisto/0.44.0

#    for R1file in `ls raw/*/$sample*_R1_001.fastq.gz`
#    do
#        flowcell=`echo $R1file | perl -ne 'print "$1" if /raw\/(.+)\//'`
#        lane=`echo $R1file | perl -ne 'print "$1" if /_(L\d\d\d)_R\d/'`
#    
#        R2file=`echo $R1file | perl -ne '($_ =~ s/_R1_/_R2_/); print $_'`
#        # Trim
#        java -XX:ParallelGCThreads=1 -Xmx2G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
#            -threads 6 \
#            -phred33 \
#            $R1file \
#            $R2file \
#            output/$sample/$flowcell.$lane.R1.fastq.gz \
#            /dev/null \
#            output/$sample/$flowcell.$lane.R2.fastq.gz \
#            /dev/null \
#            ILLUMINACLIP:ref/nebnext.fa:2:30:15:8:true \
#            TRAILING:30 \
#            MINLEN:32
#    done 
    
    # Quantify
    filelist=""
    for R1file in output/$sample/*.R1.fastq.gz
    do
        R2file=`echo $R1file | perl -ne '($_ =~ s/.R1.fastq.gz/.R2.fastq.gz/); print $_'`
        filelist="$filelist $R1file $R2file"
    done
    
    kallisto quant -i ref/Homo_sapiens.GRCh38.all.fa.idx -t 6 -o output/$sample --rf-stranded $filelist
    
    # Remove temp
    #rm output/$sample*.fastq.gz
EOF

    sbatch --export=sample=$sample --job-name=$sample.kallisto -o $script.stdout -e $script.stderr -D `pwd` $script 
    
done <unique_sample_list 

