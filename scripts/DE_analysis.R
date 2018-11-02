library(ef.utils)
library(RUVSeq)
library(reshape2)
library(dplyr)
library(ggplot2)

samples = read.table("input/samplesheet.txt", header=TRUE, sep="\t", fileEncoding="UTF-8")

libraries = Sys.glob("output/*")
libraries_label = gsub("output/(P..._V..).*", "\\1", libraries)

samples$Label = paste0(samples$Patient, "_", samples$Visite)
samples$Library = libraries[match(samples$Label, libraries_label)]
samples$kallisto = file.path(samples$Library, "abundance.tsv")
samples$TreatmentStatus = factor(ifelse(samples$Visite=="V00", "Pre", as.character(samples$Group)), levels=c("Pre", "5-ARI", "Diète"))
counts = read_identical(samples$kallisto, 1:3, 4, file.labels=samples$Label)

counts_meta = counts[,1:3]
counts_matrix = as.matrix(counts[,-(1:3)])

#filter_5 = apply(counts_matrix >=5, 1, sum) > 5
filter_20 = apply(counts_matrix >=20, 1, sum) > 20

counts_subset = counts[filter_20,]
counts_meta = counts_subset[,1:3]
counts_matrix = as.matrix(counts_subset[,-(1:3)])

ercc_genes = grepl("ERCC", counts_meta$target_id)

rownames(counts_matrix) = counts_meta$target_id
colnames(counts_matrix) = gsub(".est_counts", "", colnames(counts_matrix))
rownames(samples) = samples$Label


ercc_concentrations = read.table("input/ERCCannot.txt", sep="\t", header=TRUE)
ercc_matrix = counts_matrix[ercc_genes,]

expression_set <- newSeqExpressionSet(counts_matrix, phenoData = samples)
expression_set <- betweenLaneNormalization(expression_set, which="upper")
expression_set <- RUVg(expression_set, rownames(counts_matrix)[ercc_genes], k=1)

plotRLE(expression_set)
plotPCA(expression_set)

log_count = log2(assayData(expression_set)$normalizedCounts + 1)
log_count_pca = prcomp(t(log_count), center=FALSE, scale.=FALSE)
pca_df = as.data.frame(log_count_pca$x)
pca_df$Label = rownames(pca_df)
pca_df = dplyr::left_join(pca_df, samples, by=c(Label="Label"))
ggplot(pca_df, mapping=aes(x=PC1, y=PC2, color=Group, label=Patient)) + geom_point()

expression_subset = expression_set[,pData(expression_set)$Group!="G7"]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$TreatmentStatus = factor(expression_subset_meta$TreatmentStatus, levels=c("Pre", "5-ARI", "Diète"))
expression_subset_meta$Patient = factor(expression_subset_meta$Patient)
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta, ~ Patient + TreatmentStatus)
dds <- DESeq(dds)

