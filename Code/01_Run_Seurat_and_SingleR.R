library(rstudioapi)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(celldex)
library(SingleR)

preprocessingSCRNAseq <- function(mat, project, min.cells = 5, mt.prefix = "^mt-"){
  cat(c(">>> Project: ", project, "\n"))
  cat(c("> Input matrix: ", dim(mat), "\n"))
  cat(c("> Sum up duplicated genes\n"))
  IDfreqs <- table(rownames(mat))
  IDfreqs <- IDfreqs[IDfreqs > 1]
  if (length(IDfreqs) > 0){
    for (i in 1:length(IDfreqs)) {
      index <- which(rownames(mat) == names(IDfreqs[i]))
      duplicates <- as.matrix(mat[index, ])
      sums <- apply(duplicates, 2, sum, na.rm = TRUE)
      mat[index[1],] <- sums
      for (j in length(index):2) {
        mat <- mat[-index[j], ]
      }
    }
  }
  cat(c("matrix", dim(mat), "\n"))
  
  cat(c("> Excluding noninformative cells and features\n"))
  keep_cell <- colSums(mat > 0) > 0
  keep_feature <- rowSums(mat > 0) > 0
  mat <- mat[keep_feature, keep_cell]
  cat(c("matrix", dim(mat), "\n"))
  
  cat(c("> Check cell quality\n"))
  cell.qc <- perCellQCMetrics(mat, subsets = list(Mito = grepl(mt.prefix, rownames(mat))))
  cell.qc.res <- quickPerCellQC(cell.qc, percent_subsets = c("subsets_Mito_percent"), )
  QCMetrics <- cbind(data.frame(cell.qc), data.frame(cell.qc.res))
  ### Cell doublets or multiplets may exhibit an aberrantly high gene count
  QCMetrics$high_n_features <- isOutlier(QCMetrics$detected, nmads = 3, type = "higher", log = TRUE)
  QCMetrics <- QCMetrics %>% mutate(discard = low_lib_size | low_n_features | high_subsets_Mito_percent | high_n_features)
  rownames(QCMetrics) <- colnames(mat)
  mat.final <- mat[,!QCMetrics$discard]
  cat(c("matrix", dim(mat.final), "\n"))
  
  obj <- CreateSeuratObject(counts = mat.final, project = project, min.cells = min.cells)
  cat(c("> Create Seurat object ( min.cells =", min.cells, ")\n", dim(obj), "\n"))
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(list(SeuratObject = obj, QCMetrics = QCMetrics))
}

SeuratObjectsIntegration <- function(aggr.list){
  anchors <- FindIntegrationAnchors(object.list = aggr.list, normalization.method = "LogNormalize", 
                                    reduction = "cca", dims = 1:30)
  integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated)

  integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
  pct <- integrated@reductions[["pca"]]@stdev / sum(integrated@reductions[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)

  integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:pcs)
  integrated <- FindClusters(integrated, resolution = 0.5)

  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:pcs, n.neighbors = 15, metric = "euclidean", 
                        min.dist = 0.3)
  integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:pcs, perplexity = 30, max_iter = 3000)

  return(list(anchors = anchors, integrated = integrated, n_PCs = pcs))
}

script.folder <- paste(head(strsplit(getSourceEditorContext()$path, split="/")[[1]], -1), 
                       collapse = "/")
setwd(script.folder)
data.folder <- "../Gene_Barcode_Matrices/Int-CD45-scRNAseq_10X"
outputs.folder <- "../Outputs"

conditions <- c("WT", "KO")
ppSCRNAseq = list()
for (project.name in conditions){
  matrix.folder <- paste0(project.name, "_count")
  matrix.dir <- file.path(data.folder, matrix.folder, "filtered_gene_bc_matrices/mm10")
  matrix <- readMM(file = file.path(matrix.dir, "matrix.mtx"))
  feature.names <- read.delim(file.path(matrix.dir, "genes.tsv"), header = FALSE, stringsAsFactors = FALSE)
  barcode.names <- read.delim(file.path(matrix.dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
  colnames(matrix) <- paste(project.name, barcode.names$V1, sep = "-") %>% gsub("-", ".", .)
  rownames(matrix) <- feature.names$V2

  obj.res <- preprocessingSCRNAseq(matrix, project.name, min.cells = 2)
  saveRDS(obj.res$SeuratObject, file = file.path(outputs.folder, paste0(project.name, "_SeuratObject.rds")))
  saveRDS(obj.res$QCMetrics, file = file.path(outputs.folder, paste0(project.name, "_QCMetrics.rds")))
  ppSCRNAseq[[project.name]] <- obj.res$SeuratObject
}

comb.res <- SeuratObjectsIntegration(ppSCRNAseq)
saveRDS(comb.res$anchors, file = file.path(outputs.folder, "Seurat_integration_anchors.rds"))
saveRDS(comb.res$integrated, file = file.path(outputs.folder, "Seurat_integration_object.rds"))

integrated.RNA.counts.data <- GetAssayData(object = comb.res$integrated, slot = "counts", assay = "RNA")
file.name <- file.path(outputs.folder, "Seurat_integration_RNA_counts.txt")
write.table(integrated.RNA.counts.data, file = file.name, sep = "\t", row.names = TRUE, quote = FALSE)

PCA.cell.embeddings <- comb.res$integrated@reductions[["pca"]]@cell.embeddings
file.name <- file.path(outputs.folder, "Seurat_integration_PCA_cell_embeddings.txt")
write.table(PCA.cell.embeddings, file = file.name, sep = "\t", row.names = TRUE, quote = FALSE)

cluster.label <- as.data.frame(comb.res$integrated@active.ident)
colnames(cluster.label)[1] <- "Cluster"
add1 = function(x){
  as.integer(as.numeric(x)+1)
}
cluster.label.new = data.frame(apply(cluster.label, 2, add1), row.names = rownames(cluster.label))
file.name <- file.path(outputs.folder, "Seurat_integration_SNN_clusters.txt")
write.table(cluster.label.new, file = file.name, sep = "\t", row.names = TRUE, quote = FALSE)



# run SingleR
ref.se <- ImmGenData()
test <- as.SingleCellExperiment(comb.res$integrated, assay="RNA")
fine.single.ImmGen <- SingleR(test=test, ref=ref.se, labels=ref.se$label.fine, method="single")
labels <- as.data.frame(fine.single.ImmGen$labels, row.names=fine.single.ImmGen@rownames)
colnames(labels)[1] <- "labels"
scores <- as.data.frame(fine.single.ImmGen$scores, row.names=fine.single.ImmGen@rownames)
file.name <- file.path(outputs.folder, "SingleR_ImmGen_fine_single_labels.txt")
write.table(labels, file = file.name, sep = "\t", row.names = TRUE, quote = FALSE)
file.name <- file.path(outputs.folder, "SingleR_ImmGen_fine_single_scores.txt")
write.table(scores, file = file.name, sep = "\t", row.names = TRUE, quote = FALSE)
