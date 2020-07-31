library(RUVSeq)
library(reshape2)
library(dplyr)
library(ggplot2)
library(DESeq2)

# Determine the filename and labels for the kallisto counts for the various libraries.
libraries = Sys.glob("output/*")
libraries_label = gsub("output/(P..._V..).*", "\\1", libraries)

# Load the sample information.
samples = read.table("input/samplesheet.txt", header=TRUE, sep="\t", fileEncoding="UTF-8")
samples$Label = paste0(samples$Patient, "_", samples$Visite)
samples$Library = libraries[match(samples$Label, libraries_label)]
samples$kallisto = file.path(samples$Library, "abundance.tsv")
samples$Group = gsub("-", "", samples$Group)
samples$Group = gsub("è", "e", samples$Group)
samples$TreatmentStatus = factor(ifelse(samples$Visite=="V00", "Pre", as.character(samples$Group)), levels=c("Pre", "5ARI", "Diete"))
samples$BothConditions = samples$Patient %in% (intersect(samples %>% filter(Visite=="V00") %>% pull(Patient), samples %>% filter(Visite=="V06") %>% pull(Patient)))

# Load information about the ERCC mix added to each sample.
ercc_mix = read.table("input/ERCCmix.txt", header=TRUE, sep="\t")
samples = left_join(samples, ercc_mix, c(Label="Patient"))

# Copying read_identical from ef.utils since the package is currently impossible to load on graham due to broken dependancies.
read_identical <- function(file.names, header.columns, data.columns, file.labels=basename(file.names), sep="\t", header=TRUE, ...) {
    results=NULL
    for(i in 1:length(file.names)) {
        file.name = file.names[i]
        file.label = file.labels[i]

        file.data = read.table(file.name, sep=sep, header=header, stringsAsFactors=FALSE, ...)
        if(is.null(results)) {
            results = file.data[, header.columns, drop=FALSE]
        }

        colnames(file.data) <- paste(file.label, colnames(file.data), sep=".")

        results = cbind(results, file.data[,data.columns, drop=FALSE])
    }

    return(results)
}

# Read count information.
counts = read_identical(samples$kallisto, 1:3, 4, file.labels=samples$Label)

counts_meta = counts[,1:3]
counts_matrix = as.matrix(counts[,-(1:3)])

# Annotate transcripts.
library("biomaRt")
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations = getBM(attributes=c("ensembl_transcript_id_version", "external_gene_name", "description", "hgnc_symbol", "gene_biotype", "transcript_biotype", "transcript_tsl"),
                    filters="ensembl_transcript_id_version",
                    values=counts_meta[,1],
                    mart=ensembl)
annot_full = annotations[match(counts_meta[,1], annotations$ensembl_transcript_id_version),]
annot_full = cbind(counts_meta, annot_full[,-1])


#filter_5 = apply(counts_matrix >=5, 1, sum) > 5
filter_20 = apply(counts_matrix >=20, 1, sum) > 20

annot_full = annot_full[filter_20,]
counts_matrix = counts_matrix[filter_20,]


rownames(counts_matrix) = annot_full$target_id
colnames(counts_matrix) = gsub(".est_counts", "", colnames(counts_matrix))
rownames(samples) = samples$Label


# Load ERCC, and keep only those who should be at the same concentration in both mixes.
ercc_titration <- function(expression_data, ercc_genes, samples, label) {
    ercc_concentrations = read.table("input/ERCCannot.txt", sep="\t", header=TRUE)
    ercc_genes = annot_full$target_id %in% (ercc_concentrations %>% filter(log2.Mix.1.Mix.2.==0) %>% pull(ERCC.ID) )
    
    ercc_matrix = expression_data[ercc_genes,]
    
    ercc_df = melt(log2(ercc_matrix + 1), varnames=c("ERCC", "Patient"), value.name="Count")
    ercc_df = left_join(ercc_df, ercc_concentrations, c(ERCC="ERCC.ID"))
    ercc_df = left_join(ercc_df, samples, c(Patient="Label"))
    ercc_df$EffectiveConcentration = ifelse(ercc_df$ERCCMix==1, ercc_df$concentration.in.Mix.1..attomoles.ul., ercc_df$concentration.in.Mix.2..attomoles.ul.)
    
    ercc_mean_df = ercc_df %>% group_by(EffectiveConcentration, Patient) %>% summarize(Mean=mean(Count))
    ercc_mean_df = left_join(ercc_mean_df, samples, c(Patient="Label"))
    
    ggplot(ercc_df, aes(x=log2(EffectiveConcentration), color=Visite)) + 
        geom_point(mapping=aes(y=Count)) + 
        geom_line(ercc_mean_df, mapping=aes(y=Mean, group=Patient), alpha=0.1) +
    #    guides(color=FALSE) +
        facet_grid(~Visite)
    ggsave(paste0("output/ERCC titration ", label, ".pdf"))
}

ercc_titration(counts_matrix, ercc_genes, samples, "raw")

#ercc_concentrations = read.table("input/ERCCannot.txt", sep="\t", header=TRUE)
#ercc_genes = annot_full$target_id %in% (ercc_concentrations %>% filter(log2.Mix.1.Mix.2.==0) %>% pull(ERCC.ID) )
#
#ercc_matrix = counts_matrix[ercc_genes,]
#
#ercc_df = melt(log2(ercc_matrix + 1), varnames=c("ERCC", "Patient"), value.name="Count")
#ercc_df = left_join(ercc_df, ercc_concentrations, c(ERCC="ERCC.ID"))
#ercc_df = left_join(ercc_df, samples, c(Patient="Label"))
#ercc_df$EffectiveConcentration = ifelse(ercc_df$ERCCMix==1, ercc_df$concentration.in.Mix.1..attomoles.ul., ercc_df$concentration.in.Mix.2..attomoles.ul.)
#
#ercc_mean_df = ercc_df %>% group_by(EffectiveConcentration, Patient) %>% summarize(Mean=mean(Count))
#ercc_mean_df = left_join(ercc_mean_df, samples, c(Patient="Label"))
#
#ggplot(ercc_df, aes(x=log2(EffectiveConcentration), color=Visite)) + 
#    geom_point(mapping=aes(y=Count)) + 
#    geom_line(ercc_mean_df, mapping=aes(y=Mean, group=Patient), alpha=0.1) +
##    guides(color=FALSE) +
#    facet_grid(~Visite)
#ggsave("output/ERCC titration raw.pdf")

    
expression_set <- newSeqExpressionSet(counts_matrix, phenoData = samples)
expression_set <- betweenLaneNormalization(expression_set, which="upper")
expression_set <- RUVg(expression_set, rownames(counts_matrix)[ercc_genes], k=1)

ercc_titration(normCounts(expression_set), ercc_genes, samples, "normed")

normed_counts = cbind(normCounts(expression_set))
normed_samples = rbind(pData(expression_set))

write.table(counts_matrix, sep="\t", file="output/RawCounts.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)
write.table(normed_counts, sep="\t", file="output/NormCounts.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

# Plot MYC
myc_data = normed_counts[annot_full$external_gene_name=="MYC" & !is.na(annot_full$external_gene_name),]
myc_melt = melt(myc_data, varnames=c("Transcript", "Sample"), value.name="Count")
myc_melt = left_join(myc_melt, samples, by=c(Sample="Label"))
myc_melt$Description = factor(myc_melt$Description, levels=c("V00 Diète", "V06 Diète", "V00 5-ARI", "V06 5-ARI", "V00 G7"))
ggplot(myc_melt, aes(x=Description, y=log2(Count+1), color=Description)) + 
    geom_boxplot() + 
    facet_grid(~Transcript) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("output/MYC boxplot.pdf")
#plotRLE(expression_set)
#plotPCA(expression_set)

log_count = log2(normed_counts + 1)
log_count_pca = prcomp(t(log_count), center=FALSE, scale.=FALSE)
pca_df = as.data.frame(log_count_pca$x)
pca_df$Label = rownames(pca_df)
pca_df = dplyr::left_join(pca_df, samples, by=c(Label="Label"))

v00_df = pca_df %>% filter(Visite=="V00")
v06_df = pca_df %>% filter(Visite=="V06")
lines_df = left_join(v00_df, v06_df, by=c("Patient"="Patient"))


ggplot(pca_df, mapping=aes(x=PC1, y=PC2)) + 
    geom_point(mapping=aes(color=ERCCMix, label=Patient))

pca_df$PCA_Group = ifelse(pca_df$Description == "V00 G7", "G7", as.character(pca_df$TreatmentStatus))
ggplot(pca_df, mapping=aes(x=PC1, y=PC2)) +
    geom_segment(lines_df, mapping=aes(x=PC1.x, xend=PC1.y, y=PC2.x, yend=PC2.y), alpha=0.2) +
    geom_point(mapping=aes(color=PCA_Group, label=Patient))

ggsave("output/PCA.pdf")

# Test G7 vs rest
expression_subset = expression_set[,pData(expression_set)$TreatmentStatus =="Pre"]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$CancerStatus = factor(ifelse(expression_subset_meta$Group=="G7", "G7", "G6"), levels=c("G6", "G7"))
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta[,c("CancerStatus"), drop=FALSE], ~ CancerStatus)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(workers = 8))
save(dds, file="output/dds_group.RData")
deseq_res <- results(dds, contrast = c("CancerStatus", "G7", "G6"))
write.table(cbind(annot_full, as.data.frame(deseq_res)), file="output/CancerStatus.txt", sep="\t", col.names=TRUE, row.names=FALSE, dec=",")

# Test 5-ARI only vs rest
expression_subset = expression_set[,pData(expression_set)$Group == "5ARI" & pData(expression_set)$BothConditions]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$TreatmentStatus = factor(expression_subset_meta$TreatmentStatus, levels=c("Pre", "5ARI"))
expression_subset_meta$Patient = factor(expression_subset_meta$Patient)
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta[,c("Patient", "TreatmentStatus")], ~ Patient + TreatmentStatus)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(workers = 8))
save(dds, file="output/dds_5ARI.RData")
deseq_res <- results(dds, contrast = c("TreatmentStatus", "5ARI", "Pre"))
write.table(cbind(annot_full, as.data.frame(deseq_res)), file="output/5-ARI.txt", sep="\t", col.names=TRUE, row.names=FALSE, dec=",")

# Test Diete only vs rest
expression_subset = expression_set[,pData(expression_set)$Group =="Diete" & pData(expression_set)$BothConditions]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$TreatmentStatus = factor(expression_subset_meta$TreatmentStatus, levels=c("Pre", "Diete"))
expression_subset_meta$Patient = factor(expression_subset_meta$Patient)
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta[,c("Patient", "TreatmentStatus")], ~ Patient + TreatmentStatus)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(workers = 8))
save(dds, file="output/dds_diete.RData")
deseq_res <- results(dds, contrast = c("TreatmentStatus", "Diete", "Pre"))
write.table(cbind(annot_full, as.data.frame(deseq_res)), file="output/Diete.txt", sep="\t", col.names=TRUE, row.names=FALSE, dec=",")

# Combine both contrasts.
expression_subset = expression_set[,pData(expression_set)$Group!="G7" & pData(expression_set)$BothConditions]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$TreatmentStatus = factor(expression_subset_meta$TreatmentStatus, levels=c("Pre", "5ARI", "Diete"))
expression_subset_meta$Patient = factor(expression_subset_meta$Patient)
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta[,c("Patient", "TreatmentStatus")], ~ Patient + TreatmentStatus)
#dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- 1
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(workers = 8))
save(dds, file="output/dds.RData")

deseq_res <- results(dds, contrast = c("TreatmentStatus", "Diete", "Pre"))
write.table(as.data.frame(deseq_res), file="output/Diete.txt", sep="\t", col.names=TRUE, row.names=FALSE)
deseq_res <- results(dds, contrast = c("TreatmentStatus", "5-ARI", "Pre"))
write.table(as.data.frame(deseq_res), file="output/5-ARI.txt", sep="\t", col.names=TRUE, row.names=FALSE)


pdf("output/Venn DE.pdf")
input_files = list(Dutasteride="output/5-ARI.txt", Diete="output/Diete.txt", G7="output/CancerStatus.txt")
datasets = lapply(input_files, read.table, sep="\t", header=TRUE, quote='"', dec=",")
de_genes = lapply(datasets, function(x) { 
    x %>% filter(padj <= 0.05 & abs(log2FoldChange) > log2(1.5) & !is.na(external_gene_name)) %>% 
        pull(external_gene_name) %>% 
        unique %>% 
        as.character })
grid.draw(venn.diagram(de_genes, file=NULL))
dev.off()


#### SCHLAP1 #####
countData = read.table("output/NormCounts.txt", sep="\t", header=TRUE)
de_cancer = read.table("output/CancerStatus.txt", header=TRUE, stringsAsFactors=FALSE)
all_sample_data = read.table("input/all_samples.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

# all(rownames(countData)==de_cancer$target_id)
SCHLAP1 = which(de_cancer$external_gene_name=="SCHLAP1")
SCHLAP1_count = log2(as.numeric(countData[SCHLAP1[1],])+1)
sample_group = samples$TreatmentStatus[match(colnames(countData), samples$Label)]
gleason_score = all_sample_data$BioGleason[match(colnames(countData), all_sample_data$ID)]

schlap_df = data.frame(Expression=SCHLAP1_count,
                       Group=sample_group,
                       Gleason=as.factor(gleason_score))
                       
ggplot(schlap_df, aes(x=Gleason, y=Expression)) + geom_point()                       