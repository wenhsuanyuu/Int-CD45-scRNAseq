library(rstudioapi)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(Seurat)

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
data.folder.D6KO <- "../Gene_Barcode_Matrices/Int-CD45-scRNAseq_10X"
data.folder.Xu <- "../Gene_Barcode_Matrices/GSE124880_Xu"
outputs.folder <- "../Outputs"

### Xu et al. 2019 (GEO accession:GSE124880)
matrix <- readMM(file=file.path(data.folder.Xu, "GSE124880_PP_LP_mm10_count_matrix.mtx.gz"))
fname <- file.path(data.folder.Xu, "GSE124880_PP_LP_mm10_count_gene.tsv.gz")
feature.names <- read.delim(fname, header = FALSE, stringsAsFactors = FALSE)
fname <- file.path(data.folder.Xu, "GSE124880_PP_LP_mm10_count_barcode_v2.tsv.gz")
barcode.names <- read.delim(fname, header = FALSE, stringsAsFactors = FALSE)
colnames(matrix) <- barcode.names$V1
rownames(matrix) <- feature.names$V1

meta.data <- as.data.frame(barcode.names$V2)
colnames(meta.data) <- c("regions")
rownames(meta.data) <- barcode.names$V1

obj <- CreateSeuratObject(counts = matrix, project = "GSE124880_Xu", min.cells = 0, min.features = 0,
                          meta.data = meta.data)
ppSCRNAseq <- SplitObject(obj, split.by = "regions")

for (name in names(ppSCRNAseq)){
  if (name == "Lamina propria"){
    project.name <- "LP"
  }else{
    project.name <- "PP"
  }
  obj.res <- preprocessingSCRNAseq(as.data.frame(obj.list[[name]]@assays[["RNA"]]@counts), project.name, min.cells = 2)
  saveRDS(obj.res$SeuratObject, file = file.path(outputs.folder, paste0("GSE124880_", project.name, "_SeuratObject.rds")))
  saveRDS(obj.res$QCMetrics, file = file.path(outputs.folder, paste0("GSE124880_", project.name, "_QCMetrics.rds")))
  ppSCRNAseq[[name]] <- obj.res$SeuratObject
}

ppSCRNAseq[["Wild-type"]] = readRDS(file.path(outputs.folder, "WT_SeuratObject.rds"))
ppSCRNAseq[["DUSP6-KO"]] = readRDS(file.path(outputs.folder, "KO_SeuratObject.rds"))

comb.res <- SeuratObjectsIntegration(ppSCRNAseq)
saveRDS(comb$anchors, file = file.path(outputs.folder, "Seurat_merge_two_datasets_Anchors.rds"))
saveRDS(comb$integrated, file = file.path(outputs.folder, "Seurat_merge_two_datasets_SeuratObject.rds"))
