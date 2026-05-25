library(DeeDeeExperiment)
library(iSEE)
library(readxl)
library(S4Vectors)

#Define the following variables using the names of the files and information needed:
setwd(" ")
file_sce    <- "sce_Intri_PROVA.rds"
file_dea    <- "resultsAnnot241230.rds"
extra_columns <- c()
contrastos  <- c("Intri_C9.vs.Intri_Ctr", "Intri_C13.vs.Intri_Ctr", "Intri_C18.vs.Intri_Ctr") 
collection <- "H" 
name_file <- "dde_Intri.rds"

# In case your dataset has a "symbol" or "genenames" column, you may want to use its content as the rownames for your new DeeDeeExperiment object. 
# To do it, you just have to assign the value "TRUE" to the "convert_to" variable, and define the specific name of the column you want to use 
# in the "gene_name_col" variable.Please, be aware that in doing so, there may be some rows (genes) that will be eliminated from the final 
# object due to duplicates and empty values (NAs). 
# If you choose not to, you just have to assign the value "FALSE" to the "convert_to" variable, and in that case the object's 
# rownames will be defined by the content of the "Geneid" column. 
gene_name_col <- "GeneName"
convert_to <- TRUE 

#___________________________
sce <- readRDS(file_sce)
res_dea <- readRDS(file_dea)
gens_comuns <- intersect(rownames(res_dea), rownames(sce))

res_dea_filtrat <- res_dea[gens_comuns, ]

mostres_comunes <- intersect(colnames(sce), colnames(res_dea))

dea_list <- list()
enrich_list <- list()


for (con in contrastos) {
  #FEA
  res_fea <- paste("GSEA", collection, con, "xlsx", sep = ".")
  
  fea_up <- read_excel(res_fea, sheet = 1)
  fea_down <- read_excel(res_fea, sheet = 2)
  fea_full <- rbind(fea_up, fea_down)
  
  leading_edge_list <- replicate(nrow(fea_full), character(0), simplify = FALSE)
  
  #format 'fgsea'
  df_fea <- data.frame(
    pathway     = as.character(fea_full$ID),
    pval        = as.numeric(fea_full$pvalue),
    padj        = as.numeric(fea_full$p.adjust),
    log2err     = as.numeric(NA),
    ES          = as.numeric(fea_full$NES),
    NES         = as.numeric(fea_full$NES),
    size        = as.integer(fea_full$setSize),
    stringsAsFactors = FALSE
    )
  df_fea$leadingEdge <- leading_edge_list
  rownames(df_fea) <- df_fea$pathway
  enrich_list[[con]] <- df_fea
  
  #DEA
  df_dea <- res_dea_filtrat
  if (length(contrastos) > 1) {
    altres_cons <- contrastos[contrastos != con] #contrastos que NO mirem ARA
  
    pattern_excloure <- paste(altres_cons, collapse = "|") 
    cols_keep <- !grepl(pattern_excloure, colnames(res_dea_filtrat)) #agafa tot el que NO esta dins de pattern_excloure
  
    df_dea <- res_dea_filtrat[, cols_keep]
  }
  colnames(df_dea)[colnames(df_dea) == paste0("logFC.", con)]     <- "log2FoldChange"
  colnames(df_dea)[colnames(df_dea) == paste0("P.Value.", con)]   <- "pvalue"
  colnames(df_dea)[colnames(df_dea) == paste0("adj.P.Val.", con)] <- "padj"
  colnames(df_dea)[colnames(df_dea) == paste0("FC.", con)] <- "FoldChange"
  
  dea_list[[con]] <- as.data.frame(df_dea)
}

dde <- DeeDeeExperiment(sce = sce)

#lost cols 
cols_interes <- c("Geneid", "Symbol", "Chr", "Start", "End", "Strand", "length", "Description", extra_columns)
cols_existents <- unlist(lapply(cols_interes, function(x) {
  grep(paste0("^", x, "$"), colnames(res_dea_filtrat), ignore.case = TRUE)
}))

info_extra <- res_dea_filtrat[rownames(dde), cols_existents, drop = FALSE]
for (col in colnames(info_extra)) {
  rowData(dde)[[col]] <- info_extra[[col]]
}

col_symbol <- grep(paste0("^symbol$|", gene_name_col), colnames(rowData(dde)), ignore.case = TRUE, value = TRUE)

if (length(col_symbol) > 0 && convert_to == TRUE) {
  if (length(col_symbol) > 1) {
    col_symbol <- col_symbol[1]
  }
  actual_symbol <- rowData(dde)[[col_symbol]]
  keep <- !is.na(actual_symbol) & actual_symbol != "" & !duplicated(actual_symbol)
  
  dde <- dde[keep, ] #eliminem files de counts i rowData alhora
  nous_noms <- rowData(dde)[[col_symbol]]
  rownames(dde) <- nous_noms
  
  for (con in names(dea_list)) {
    dea_list[[con]] <- dea_list[[con]][keep, ]
    rownames(dea_list[[con]]) <- nous_noms
  }
}

dde <- addDEA(dde, dea = dea_list)
dde <- addFEA(dde, fea = enrich_list)

saveRDS(dde, file = name_file)

