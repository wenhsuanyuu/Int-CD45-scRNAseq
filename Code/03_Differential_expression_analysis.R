library(rstudioapi)
library(limma)
library(edgeR)

script.folder <- paste(head(strsplit(getSourceEditorContext()$path, split="/")[[1]], -1), 
                       collapse = "/")
setwd(script.folder)
inputs.folder <- "../Outputs/KOvWT_DEG_inputs"
results.folder <- "../Outputs/KOvWT_DEG_results"

celltypes <- c("B cell", "CD4+ T cell", "CD8+ T cell", "DC", "ILC1", "ILC2", 
               "Macrophage", "Monocyte", "NK cell", "NKT", "Tgd", "Treg")
for (celltype in celltypes){
  fname <- paste("Differential_expression_analysis", celltype, "matrix.txt", sep = "_")
  x <- read.delim(file.path(inputs.folder, fname), sep = "\t", row.names=1)
  fname <- paste("Differential_expression_analysis", celltype, "group.txt", sep = "_")
  group <- scan(file.path(inputs.folder, fname), what=character(), sep="\t")
  
  dge <- DGEList(x, group = group)
  dge <- calcNormFactors(dge)
  cdr <- scale(colMeans(x > 0))
  design <- model.matrix(~ cdr + group)
  rownames(design) <- colnames(x)
  dge <- estimateDisp(dge, design=design, robust=TRUE)
  fit <- glmQLFit(dge, design=design)
  qlf <- glmQLFTest(fit)
  
  temp <- as.data.frame(qlf$AveLogCPM)
  colnames(temp) <- c("AveLogCPM")
  rownames(temp) <- rownames(qlf$table)
  fname = paste("Differential_expression_analysis", celltype, "AveLogCPM.txt", sep = "_")
  write.table(temp, file=file.path(results.folder, fname), sep="\t", row.names=TRUE, col.names=TRUE, 
              quote = FALSE)
  
  tt <- topTags(qlf, n=Inf)
  results <- as.data.frame(tt)
  fname = paste("Differential_expression_analysis", celltype, "results.txt", sep = "_")
  write(paste("GeneName", "logFC", "logCPM", "F", "PValue", "FDR", sep="\t"), file=fname)
  write.table(results, file=file.path(results.folder, fname), append=TRUE, sep="\t", 
              row.names=TRUE, col.names=FALSE, quote = FALSE)
}

