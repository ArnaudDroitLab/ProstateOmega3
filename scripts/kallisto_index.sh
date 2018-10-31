module load kallisto/0.44.0

# Generate ERCC fasta
mkdir output
awk '{print ">" $1; print $5;}' <  input/ERCCseq.txt >  output/ERCC.fa

# Download ensembl reference transcriptome.
mkdir ref
cd ref
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..

# Generate kalisto index
kallisto index -i ref/Homo_sapiens.GRCh38.all.fa.idx ref/Homo_sapiens.GRCh38.cdna.all.fa ref/HS_rRNA.fa ref/Homo_sapiens.GRCh38.ncrna.fa  output/ERCC.fa

