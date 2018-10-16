module load dataset/Homo_sapiens_GRCh38
module load trimmomatic/0.36
module load kallisto/0.44.0

# Download ensembl reference transcriptome.
mkdir ref
cd ref
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# Generate ERCC fasta
awk '{print ">" $1; print $5;}' < input/ERCCseq.txt > output/ERCC.fa

# Generate kalisto index
kallisto index -i Homo_sapiens.GRCh38.cdna.all.fa.idx Homo_sapiens.GRCh38.cdna.all.fa ../output/ERCC.fa
cd ..

# Make a list of all unique samples.
rm tmp_file_list
curdir=`pwd`
for i in raw/180925_D00487S_0161_ACCR9DANXX/ raw/180925_D00487S_0162_BCCU7UANXX/ raw/181003_D00487S_0163_ACCU7FANXX/
do
    pushd $i
    ls * >> $curdir/tmp_file_list
    popd
done

# Make a list of all samples.
sort < tmp_file_list | perl -ne 'print "$1\n" if /(.+)_L..._R._001.fastq.gz/' | sort | uniq > unique_sample_list

for read sample; do
    mkdir -p output/$sample
    
    for R1file in `ls /u01/genomique/Demultiplexing/results/*/fastq/20180703-BILS001/$sample*_R1_001.fastq.gz`
    do
        flowcell=`echo $R1file | perl -ne 'print "$1" if /results\/(.+)\/fastq/'`
        lane=`echo $R1file | perl -ne 'print "$1" if /_(L\d\d\d)_R\d/'`
    
        R2file=`echo $R1file | perl -ne '($_ =~ s/_R1_/_R2_/); print $_'`
        # Trim
        java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
            -threads 6 \
            -phred33 \
            $R1file \
            $R2file \
            output/$sample/$flowcell.$lane.R1.fastq.gz \
            /dev/null \
            output/$sample/$flowcell.$lane.R2.fastq.gz \
            /dev/null \
            ILLUMINACLIP:ref/nebnext.fa:2:30:15:8:true \
            TRAILING:30 \
            MINLEN:32
    done 
    
    # Quantify
    filelist=""
    for R1file in output/$sample/*.R1.fastq.gz
    do
        R2file=`echo $R1file | perl -ne '($_ =~ s/.R1.fastq.gz/.R2.fastq.gz/); print $_'`
        filelist="$filelist $R1file $R2file"
    done
    
    kallisto quant -i ref/Homo_sapiens.GRCh38.cdna.all.fa.idx -o output/$sample --rf-stranded $filelist
    
    # Remove temp
    # rm output/$sample*.fastq.gz
    
done <unique_sample_list 




JOB_NAME=trimmomatic.DEX_nonMamm_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_nonMamm_1_1.07127a1c5c6871e58320f42ab51a3420.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_nonMamm_1_1.07127a1c5c6871e58320f42ab51a3420.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.log
trimmomatic.DEX_nonMamm_1_1.07127a1c5c6871e58320f42ab51a3420.mugqic.done
)