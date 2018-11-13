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
    script=output/$sample/align_star.sh
    cat <<'EOF' > $script
#!/bin/bash
#SBATCH --account=def-stbil30
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

module load star/2.6.1a
    
    # Quantify
    filelist=""
    R1files=""
    R2files=""
    for R1file in output/$sample/*.R1.fastq.gz
    do
        R1files="$R1files,$R1file"
        
        R2file=`echo $R1file | perl -ne '($_ =~ s/.R1.fastq.gz/.R2.fastq.gz/); print $_'`
        R2files="$R2files,$R2file"
    done
    
    R1files=`echo $R1files | sed -e 's/^,//'`
    R2files=`echo $R2files | sed -e 's/^,//'`
    filelist="$R1files $R2files"
    
    STAR --runMode alignReads \
    --genomeDir ref/STAR_index \
    --readFilesIn $filelist \
    --runThreadN 14 \
    --readFilesCommand zcat \
    --outStd Log \
    --outSAMunmapped Within \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix output/$sample/STAR \
    --limitGenomeGenerateRAM 100000000000 \
    --limitIObufferSize 4000000000

    
    # Remove temp
    #rm output/$sample*.fastq.gz
EOF

    sbatch --export=sample=$sample --job-name=$sample.STAR -o $script.stdout -e $script.stderr -D `pwd` $script 
    
done <unique_sample_list 




