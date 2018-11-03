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
samples$TreatmentStatus = factor(ifelse(samples$Visite=="V00", "Pre", as.character(samples$Group)), levels=c("Pre", "5-ARI", "Diète"))
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

#filter_5 = apply(counts_matrix >=5, 1, sum) > 5
filter_20 = apply(counts_matrix >=20, 1, sum) > 20

counts_subset = counts[filter_20,]
counts_meta = counts_subset[,1:3]
counts_matrix = as.matrix(counts_subset[,-(1:3)])


rownames(counts_matrix) = counts_meta$target_id
colnames(counts_matrix) = gsub(".est_counts", "", colnames(counts_matrix))
rownames(samples) = samples$Label


ercc_concentrations = read.table("input/ERCCannot.txt", sep="\t", header=TRUE)
ercc_genes = counts_meta$target_id %in% (ercc_concentrations %>% filter(log2.Mix.1.Mix.2.==0) %>% pull(ERCC.ID) )

ercc_matrix = counts_matrix[ercc_genes,]

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

# We'll normalize in two batches: one for mix1 data, and the other for mix 2 data.
# That means splitting V00 from V06 data.
    
expression_set <- newSeqExpressionSet(counts_matrix, phenoData = samples)
expression_set <- betweenLaneNormalization(expression_set, which="upper")
expression_set <- RUVg(expression_set, rownames(counts_matrix)[ercc_genes], k=1)

normed_counts = cbind(normCounts(expression_set))
normed_samples = rbind(pData(expression_set))


#plotRLE(expression_set)
#plotPCA(expression_set)

log_count = log2(normed_counts + 1)
log_count_pca = prcomp(t(log_count), center=FALSE, scale.=FALSE)
pca_df = as.data.frame(log_count_pca$x)
pca_df$Label = rownames(pca_df)
pca_df = dplyr::left_join(pca_df, samples, by=c(Label="Label"))
ggplot(pca_df, mapping=aes(x=PC1, y=PC2, color=ERCCMix, label=Patient)) + geom_point()

expression_subset = expression_set[,pData(expression_set)$Group!="G7" & pData(expression_set)$BothConditions & pData(expression_set)$Patient %in% c("P031", "P045", "P049", "P050", "P053")]
expression_subset_meta = pData(expression_subset)
expression_subset_meta$TreatmentStatus = factor(expression_subset_meta$TreatmentStatus, levels=c("Pre", "5-ARI", "Diète"))
expression_subset_meta$Patient = factor(expression_subset_meta$Patient)
dds <- DESeqDataSetFromMatrix(normCounts(expression_subset), expression_subset_meta, ~ Patient + TreatmentStatus)
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(workers = cores))
save(dds, file="output/dds.RData")

deseq_res <- results(dds, contrast = c("TreatmentStatus", "Diète", "Pre"))

deseq_res <- results(dds, contrast = c("TreatmentStatus", "5-ARI", "Pre"))