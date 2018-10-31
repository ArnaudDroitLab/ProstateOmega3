samples = read.table("input/samplesheet.txt", header=TRUE, sep="\t")

libraries = Sys.glob("output/*")
libraries_label = gsub("output/(P..._V..).*", "\\1", libraries)

samples$label = paste0(samples$Patient, "_", samples$Visite)
samples$libraries = libraries[match(samples$label, libraries_label)]
samples$kallisto = file.path(samples$libraries, "abundance.tsv")