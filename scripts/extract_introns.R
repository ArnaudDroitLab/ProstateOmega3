library(GenomicRanges)
library(dplyr)

ensembl_full_gtf = rtracklayer::import("ref/Homo_sapiens.GRCh38.94.chr.gtf")

ensembl_gene_gtf = ensembl_full_gtf[ensembl_full_gtf$type=="gene"]
ensembl_exon_gtf = ensembl_full_gtf[ensembl_full_gtf$type=="exon"]

introns = GenomicRanges::setdiff(ensembl_gene_gtf, ensembl_exon_gtf)
intron_gene_match = findOverlaps(introns, ensembl_gene_gtf)
introns$source[queryHits(intron_gene_match)]= as.character(ensembl_gene_gtf$source[subjectHits(intron_gene_match)])


introns$gene_id[queryHits(intron_gene_match)]= paste0(as.character(ensembl_gene_gtf$gene_id[subjectHits(intron_gene_match)]), "_intron")
introns$transcript_id = paste0(introns$gene_id, "_transcript")
#introns$gene_id[queryHits(intron_gene_match)]= paste0(as.character(ensembl_gene_gtf$gene_id[subjectHits(intron_gene_match)]), "_intron_", queryHits(intron_gene_match))
introns$type="exon"

# Set exon number.
i=1
exon_number=NA
gene_ids = introns$gene_id
while(any(is.na(exon_number))) {
    current_exons = duplicated(gene_ids)
    exon_number[!current_exons] = i
    gene_ids[!current_exons] = "DUPLICATE"
    
    i = i + 1
}
# The very first entry will always be overwritten since 
# it is the first DUPLICATE dummy entry. Being the first entry, it is also
# always a exon number 1.
exon_number[1] = 1
introns$exon_number = exon_number

# Add gene and transcript entries.
intron_transcripts = as.data.frame(introns) %>% group_by(gene_id) %>% summarize(seqnames=seqnames[1], start=min(start), end=max(end), strand=strand[1], source=source[1])
intron_transcripts$start =as.integer(intron_transcripts$start)
intron_transcripts$end =as.integer(intron_transcripts$end)
intron_transcripts$exon_number = NA
intron_transcripts$type ="transcript"

intron_genes = intron_transcripts
intron_genes$type = "gene"

intron_transcripts$transcript_id = paste0(intron_transcripts$gene_id, "_transcript")
intron_genes$transcript_id = NA

intron_transcripts = intron_transcripts[,c("seqnames", "start", "end", "strand", colnames(mcols(introns)))]
intron_genes = intron_genes[,c("seqnames", "start", "end", "strand", colnames(mcols(introns)))]

all_intron_features = c(introns, GRanges(intron_transcripts), GRanges(intron_genes))
feature_type_order = ifelse(all_intron_features$type=="gene", 1, ifelse(all_intron_features$type=="transcript", 2, all_intron_features$exon_number + 2))
feature_order = order(all_intron_features$gene_id, feature_type_order)

all_intron_features = all_intron_features[feature_order]

ensembl_colsubset = ensembl_full_gtf
mcols(ensembl_colsubset) = mcols(ensembl_colsubset)[colnames(mcols(introns))]

all_features=c(ensembl_colsubset, all_intron_features)

rtracklayer::export(all_features, con="output/introns_combined.gtf", format="gtf")