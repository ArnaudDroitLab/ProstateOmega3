library(ef.utils)
library(RUVSeq)

samples = read.table("input/samplesheet.txt", header=TRUE, sep="\t")

libraries = Sys.glob("output/*")
libraries_label = gsub("output/(P..._V..).*", "\\1", libraries)

samples$label = paste0(samples$Patient, "_", samples$Visite)
samples$libraries = libraries[match(samples$label, libraries_label)]
samples$kallisto = file.path(samples$libraries, "abundance.tsv")

counts = read_identical(samples$kallisto, 1:3, 4, file.labels=samples$label)

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
rownames(samples) = samples$label


ercc_concentrations = read.table("input/ERCCannot.txt", sep="\t", header=TRUE)
ercc_matrix = counts_matrix[ercc_genes,]

expression_set <- newSeqExpressionSet(counts_matrix, phenoData = samples)
expression_set <- betweenLaneNormalization(expression_set, which="upper")
expression_set <- RUVg(expression_set, rownames(counts_matrix)[ercc_genes], k=1)

plotRLE(expression_set_RUVg)
plotPCA(expression_set_RUVg)