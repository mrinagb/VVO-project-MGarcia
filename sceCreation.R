library(SingleCellExperiment)
library(readxl)
library(DESeq2)

#Define the following variables using the names of the files and information needed:
dataDir <- setwd("")
# load(file.path(dataDir,"deseq2.dds.RData")) 
# counts_matrix <- counts(dds)
# counts_matrix <- readRDS("rawCounts_allbatches.rds")
counts <- read.table(file = file.path(dataDir, "gene_count.xls"), sep = "\t", header = T, dec = ".", stringsAsFactors = F)
rownames(counts) <- counts$gene_id
counts_matrix <- counts[,2:31]

sample_metadata <- read.table("Targets.txt", header = TRUE, sep ="\t", check.names = FALSE)
gene_metadata_raw <- as.data.frame(read_excel("resultsDEA_20241230SV_SYM.xlsx"))
extra_columns <- "GeneName"
name_file <- "sce_Intri_PROVA.rds"

#duplicates filter
gene_metadata_raw <- gene_metadata_raw[!duplicated(gene_metadata_raw$Geneid), ]
rownames(gene_metadata_raw) <- gene_metadata_raw$Geneid

rownames(sample_metadata) <- sample_metadata$sampleName #sampleName

wishcols <- c("Geneid", "Symbol", "Chr", "Start", "End", "Strand", "Description", "length", extra_columns)
cols_trobades <- unlist(lapply(wishcols, function(x) {
  grep(paste0("^", x, "$"), colnames(gene_metadata_raw), ignore.case = TRUE)
}))

gene_metadata_clean <- gene_metadata_raw[, cols_trobades]

#common genes
gens_comuns <- intersect(rownames(counts_matrix), rownames(gene_metadata_clean))
gene_mtdt_final <- gene_metadata_clean[gens_comuns, ]

#common samples  
mostres_comuns <- intersect(colnames(counts_matrix), rownames(sample_metadata))
counts_final <- counts_matrix[gens_comuns, mostres_comuns] #+gens_comuns!!!!
sample_mtdt_final <- sample_metadata[mostres_comuns, ]

stopifnot(nrow(counts_final) == nrow(gene_mtdt_final))
stopifnot(identical(rownames(counts_final), rownames(gene_mtdt_final)))
stopifnot(ncol(counts_final) == nrow(sample_mtdt_final))

sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(counts_final)),
  colData = sample_mtdt_final,
  rowData = gene_mtdt_final
)

saveRDS(sce, file = name_file)
