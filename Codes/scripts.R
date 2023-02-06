# Wilms Reanalysis
# work.path <- "/share/genomics/cwx/External/xlu/wilms";setwd(work.path)
work.path <- "i:/genomicdata/External/Xlu/Wilms";setwd(work.path)
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
code.path <- file.path(work.path, "Codes")
fig.path <- file.path(work.path, "Figures")
report.path <- file.path(work.path, "Reports")

invisible(lapply(ls()[grep("path", ls())], function(x){
  if (!dir.exists(get(x))) dir.create(get(x))
}))


library(Matrix)
library(Seurat)
library(magrittr)
library(SingleR)
library(openxlsx)
library(harmony)

# ##
reference <- readRDS("/share/genomics/cwx/Public/reference.rds")
# hpca.se <- reference$HumanPrimaryCellAtlasData
# blueprint <- reference$BlueprintEncodeData
# immune.se <- reference$DatabaseImmuneCellExpressionData

LM22 <- read.table("LM22", sep = "\t", row.names = 1, header = T)
LM22 <- LM22[rownames(LM22)%in%rownames(seu), ]
LM22.signature <- lapply(colnames(LM22), function(x){
  x = setNames(object = LM22[[x]], nm = rownames(LM22))
  names(sort(x, decreasing = T))[1:50]
})
names(LM22.signature) <- colnames(LM22)


# Data Preparation --------------------------------------------------------

## expression --------------------------------------------------------
barcodes <- read.table(file = file.path(data.path, "tableOfCounts_colLabels.tsv"),
                       header = T, sep = "\t")
genes <- read.table(file = "data/tableOfCounts_rowLabels.tsv",
                    header = T, sep = "\t")
mtx <- readMM(file = "data/tableOfCounts.mtx")
mtx <- as(mtx,"data/dgCMatrix")
colnames(mtx) <- barcodes$DropletID; rownames(mtx) <- genes$Symbol
saveRDS(mtx, "raw.expr.rds")
mtx <- mtx[rowSums(mtx)>0,]
sum(duplicated(rownames(mtx)));rownames(mtx)[duplicated(rownames(mtx))]
# mtx <- mtx[rownames(mtx)[!duplicated(rownames(mtx))],rownames(metadata)]

dim(mtx)
summary(rowSums(mtx))
summary(rowSums(mtx>0))
summary(colSums(mtx))
summary(colSums(mtx>0))

## cell info --------------------------------------------------------
sheets.name <- getSheetNames(file.path(data.path, "aat1699-Young-TablesS1-S12-revision2.xlsx"))
sheets <- lapply(sheets.name, function(x){
  read.xlsx(file.path(data.path, "aat1699-Young-TablesS1-S12-revision2.xlsx"), sheet = x)
})
sheets <- lapply(sheets, as.data.frame)
names(sheets) <- sheets.name

study_info <- sheets$`TableS6 - Sample manifest`
colnames(study_info) <- study_info[1,];study_info <- study_info[-1,]
cell_info <- sheets$`TableS11 - Cell manifest`
colnames(cell_info) <- cell_info[1,];cell_info <- cell_info[-1,]
cell_info <- subset(cell_info, QCpass == "TRUE")
cluster_info <- sheets$`TableS2 - Cluster info`
colnames(cluster_info) <- cluster_info[1, ]; cluster_info <- cluster_info[-1, ]
rownames(cluster_info) <- cluster_info$Cluster_ID
cluster_info <- cluster_info[c("Category", "Cell_type1", 
                               "Cell_type2", "Cell_type3", "Genotype")]

metadata <- cbind(cell_info, study_info[match(cell_info$SangerID,study_info$SangerID), ])
metadata <- metadata[c("barcode","SangerID","ClusterID","Compartment","Source","Channel10X",
                       "Experiment","TissueDiseaseState","Organ","Location1","Location2",
                       "BiologicalRepNo","TechnicalRepNo","Sort","AgeInMonthsPostConception",
                       "PatientDiseaseState")]
metadata <- cbind(metadata, cluster_info[match(metadata$ClusterID, rownames(cluster_info)), ])
rownames(metadata) <- cell_info$DropletID
saveRDS(metadata, file.path(res.path, "metadata.rds"))


# Clustering and Annotation -----------------------------------------------


## Clustering --------------------------------------------------------------

mtx <- readRDS(file.path(data.path, "raw.expr.rds"))
all(rownames(metadata)%in%colnames(mtx))
mtx <- mtx[,rownames(metadata)]

## Preprocess
seu <- CreateSeuratObject(counts = mtx, min.cells = 3)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- AddMetaData(seu, metadata = metadata)
seu <- subset(seu, PatientDiseaseState == "Wilms")
Idents(seu) <- seu$Experiment
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
# seu <- subset(seu, nFeature_RNA >= 200 & nFeature_RNA <=6000 & percent.mt < 30)
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- RunPCA(seu)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:15, reduction.name = "umap_pca")
DimPlot(seu, reduction = "pca", split.by = "Channel10X", ncol = 4)
DimPlot(seu, reduction = "umap_pca", split.by = "Channel10X", ncol = 4)

## Remove Batch Effect
seu <- RunHarmony(seu, group.by.vars = "Channel10X")
# seu <- RunHarmony(seu, group.by.vars = "Patient", assay.use = "SCT")
ElbowPlot(seu, reduction = "harmony", ndims = 50)
dim.to.use = order(seu@reductions$harmony@stdev, decreasing = T)[1:15]
seu <- RunUMAP(seu, reduction = "harmony", dims = dim.to.use, reduction.name = "umap_harmony")

DimPlot(seu, reduction = "umap_harmony", split.by = "Channel10X", ncol = 4)
DimPlot(seu, reduction = "umap_harmony", group.by = "Channel10X")
DimPlot(seu, reduction = "umap_harmony", group.by = "Category")
DimPlot(seu, reduction = "umap_harmony", group.by = "Cell_type1")
DimPlot(seu, reduction = "umap_harmony", group.by = "Cell_type2")
DimPlot(seu, reduction = "umap_harmony", group.by = "Cell_type3")
DimPlot(seu, reduction = "umap_harmony", group.by = "seurat_clusters")

plot.marker <- c("PECAM1", "VWF", "EPCAM", "KRT19", 
                 "PTPRC", "CD3D", "CD4", "CD8A", "NKG7", 
                 "CD79A", "MS4A1", 
                 "CD14", "MRC1", "S100A8",
                 "PTGDS", "DCN", "LUM", "ACTA2", 
                 "TPSAB1", "TPSB2")
FeaturePlot(seu, features = plot.marker, reduction = "umap_harmony")
FeaturePlot(seu, features = plot.marker, reduction = "umap_pca")

seu <- FindNeighbors(seu, reduction = "harmony", dims = dim.to.use)
seu <- FindClusters(seu, resolution = 0.4)

markers <- FindAllMarkers(seu)
markers <- markers[markers$avg_log2FC>0 & markers$p_val_adj<0.05, ]
View(markers[markers$cluster == 4, ])
saveRDS(markers, file.path(res.path, "markers.rds"))

## Annotation --------------------------------------------------------------

DimPlot(seu, group.by = "seurat_clusters", reduction = "umap_harmony", label = T)
if (!file.exists(file.path(res.path, "anno.txt"))){
  write.table(data.frame(label = levels(seu$seurat_clusters), celltypes = "unknown"), 
              file.path(res.path, "anno.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}
anno <- read.table(file.path(res.path, "anno.txt"), sep = "\t", header = T)
seu$celltypes <- plyr::mapvalues(x = seu$seurat_clusters,
                                 from = anno$label,
                                 to = anno$celltypes)
DimPlot(seu, group.by = "celltypes", reduction = "umap_harmony", label = T)

ref <- reference$HumanPrimaryCellAtlasData
SingleR.anno <- list(
  "main" = as.data.frame(SingleR(test = seu@assays$RNA@data, 
                                 ref = ref, 
                                 labels = ref$label.main)),
  "fine" = as.data.frame(SingleR(test = seu@assays$RNA@data, 
                                 ref = ref, 
                                 labels = ref$label.fine))
)

seu$SingleR.main.anno <- SingleR.anno$main$pruned.labels
seu$SingleR.fine.anno <- SingleR.anno$fine$pruned.labels

saveRDS(seu, file.path(res.path, "seu.rds"))


# T cell Annotation -------------------------------------------------------

T.seu <- subset(seu, celltypes %in% c("T cell", "NK cell"))
T.seu <- FindVariableFeatures(T.seu) %>% ScaleData()
T.seu <- RunPCA(T.seu)
T.seu <- RunUMAP(T.seu, dims = 1:30, reduction.name = "umap_pca")
DimPlot(T.seu, reduction = "pca", split.by = "Channel10X", ncol = 4)
DimPlot(T.seu, reduction = "umap_pca", group.by = "Channel10X")
DimPlot(T.seu, reduction = "umap_pca")

## Remove Batch Effect
T.seu <- RunHarmony(T.seu, group.by.vars = "Channel10X", max.iter.harmony = 50)
ElbowPlot(T.seu, reduction = "harmony", ndims = 50)
dim.to.use = order(T.seu@reductions$harmony@stdev, decreasing = T)[1:30]
T.seu <- RunUMAP(T.seu, reduction = "harmony", dims = dim.to.use, reduction.name = "umap_harmony")

DimPlot(T.seu, reduction = "umap_harmony", split.by = "Channel10X", ncol = 4)
DimPlot(T.seu, reduction = "umap_harmony", group.by = "Channel10X")
DimPlot(T.seu, reduction = "umap_harmony", group.by = "Category")
DimPlot(T.seu, reduction = "umap_harmony", group.by = "Cell_type1")
DimPlot(T.seu, reduction = "umap_harmony", group.by = "Cell_type2")
DimPlot(T.seu, reduction = "umap_harmony", group.by = "Cell_type3")
DimPlot(T.seu, reduction = "umap_harmony", group.by = "seurat_clusters", label = T)

T.seu <- FindNeighbors(T.seu, reduction = "harmony", dims = dim.to.use)
T.seu <- FindClusters(T.seu, resolution = 1.2)
T.avg.expr <- AverageExpression(T.seu, group.by = "seurat_clusters")[[1]][rownames(LM22), ]
pheatmap(cor(LM22[, 4:12], T.avg.expr, method = "spearman"), display_numbers = T)

T.markers <- FindAllMarkers(T.seu)
T.markers <- T.markers[T.markers$avg_log2FC>0 & T.markers$p_val_adj<0.05, ]
View(T.markers[T.markers$cluster == 8, ])
saveRDS(T.markers, file.path(res.path, "T.markers.rds"))

## Annotation From Raw Article ---------------------------------------------


## Marker-based Annotation -------------------------------------------------

FeaturePlot(T.seu, features = c("CD3D", "CD4", "CD8A", "CD8B", "CD38",
                                "FOXP3", "IL2RA", "XCL1", "LEF1", "CCR7"),
            reduction = "umap_harmony")
DimPlot(T.seu, group.by = "seurat_clusters", reduction = "umap_harmony", label = T)
if (!file.exists(file.path(res.path, "T.anno.txt"))){
  write.table(data.frame(label = levels(T.seu$seurat_clusters), celltypes = "unknown"), 
              file.path(res.path, "T.anno.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}
T.anno <- read.table(file.path(res.path, "T.anno.txt"), sep = "\t", header = T)
T.seu$celltypes <- plyr::mapvalues(x = T.seu$seurat_clusters,
                                 from = T.anno$label,
                                 to = T.anno$celltypes)
DimPlot(T.seu, group.by = "celltypes", reduction = "umap_harmony")


## SingleR Annotation ------------------------------------------------------

T.SingleR.anno <- list(
  "main" = as.data.frame(SingleR(test = T.seu@assays$RNA@data, 
                                 ref = reference$DatabaseImmuneCellExpressionData, 
                                 labels = reference$DatabaseImmuneCellExpressionData$label.main)),
  "fine" = as.data.frame(SingleR(test = T.seu@assays$RNA@data, 
                                 ref = reference$DatabaseImmuneCellExpressionData, 
                                 labels = reference$DatabaseImmuneCellExpressionData$label.fine))
)
T.seu$SingleR.main.anno <- T.SingleR.anno$main$pruned.labels
T.seu$SingleR.fine.anno <- T.SingleR.anno$fine$pruned.labels
DimPlot(T.seu, group.by = "SingleR.main.anno", reduction = "umap_harmony")
DimPlot(T.seu, group.by = "SingleR.fine.anno", reduction = "umap_harmony")
saveRDS(T.seu, file.path(res.path, "T.seu.rds"))


# Myeloid cell ------------------------------------------------------------------

Mye.seu <- subset(seu, celltypes == "Myeloid")
Mye.seu <- FindVariableFeatures(Mye.seu) %>% ScaleData()
Mye.seu <- RunPCA(Mye.seu)
Mye.seu <- RunUMAP(Mye.seu, dims = 1:30, reduction.name = "umap_pca")
DimPlot(Mye.seu, reduction = "pca", split.by = "Channel10X", ncol = 4)
DimPlot(Mye.seu, reduction = "umap_pca", split.by = "Channel10X", ncol = 4)
DimPlot(Mye.seu, reduction = "umap_pca", group.by = "Channel10X")
DimPlot(Mye.seu, reduction = "umap_pca")

## Remove Batch Effect
Mye.seu <- RunHarmony(Mye.seu, group.by.vars = "Channel10X", max.iter.harmony = 50)
ElbowPlot(Mye.seu, reduction = "harmony", ndims = 50)
dim.to.use = order(Mye.seu@reductions$harmony@stdev, decreasing = T)[1:30]
Mye.seu <- RunUMAP(Mye.seu, reduction = "harmony", dims = dim.to.use, reduction.name = "umap_harmony")

DimPlot(Mye.seu, reduction = "umap_harmony", split.by = "Channel10X", ncol = 4)
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "Channel10X")
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "Category")
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "Cell_type1")
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "Cell_type2")
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "Cell_type3")
DimPlot(Mye.seu, reduction = "umap_harmony", group.by = "seurat_clusters")

Mye.seu <- FindNeighbors(Mye.seu, reduction = "harmony", dims = dim.to.use)
Mye.seu <- FindClusters(Mye.seu, resolution = 1.2)

FeaturePlot(Mye.seu, reduction = "umap_pca",
            features =  c("CD68", "CD86", "TLR4", "VSIR", "CD80", "CCR7", "CD163",
                          "CD14", "CD68", "CD16", "LYZ", "FCN1", "S100A8", "S100A9",
                          "CD1C", "FCER1A", "S100A8", "S100A9", "MS4A1", "CD79A"))
FeaturePlot(Mye.seu, reduction = "umap_harmony",
            features =  c("CD14", "CD16", "CD86", "LYZ", "FCN1", "VCAN", "S100A8", "S100A9",
                          "CD68", "CD86", "TLR4", "VSIR", "CCR7", "CD163", "CD1C",
                          "FCER1A", "CLEC9A", "CLEC10A", "FLT3", "CD79A", "MS4A1"))

Mye.markers <- FindAllMarkers(Mye.seu)
Mye.markers <- Mye.markers[Mye.markers$avg_log2FC>0&Mye.markers$p_val_adj<0.05, ]
View(Mye.markers[Mye.markers$cluster == 8, ])
Mye.markers$anno <- sapply(Mye.markers$gene, LM22.annot)
Mye.avg.expr <- AverageExpression(Mye.seu, group.by = "seurat_clusters")[[1]][rownames(LM22), ]
pheatmap(cor(LM22[, c(1:3, 13:22)], Mye.avg.expr, method = "spearman"), display_numbers = T)
saveRDS(Mye.markers, file.path(res.path, "Mye.markers.rds"))

DimPlot(Mye.seu, group.by = "seurat_clusters", reduction = "umap_harmony", label = T)
if (!file.exists(file.path(res.path, "Mye.anno.txt"))){
  write.table(data.frame(label = levels(Mye.seu$seurat_clusters), celltypes = "unknown"), 
              file.path(res.path, "Mye.anno.txt"),
              row.names = F, col.names = T, quote = F, sep = "\t")
}
Mye.anno <- read.table(file.path(res.path, "Mye.anno.txt"), sep = "\t", header = T)
Mye.seu$celltypes <- plyr::mapvalues(x = Mye.seu$seurat_clusters,
                                   from = Mye.anno$label,
                                   to = Mye.anno$celltypes)
DimPlot(Mye.seu, group.by = "celltypes", reduction = "umap_harmony")

ref <- reference$HumanPrimaryCellAtlasData
ref <- ref[, ref$label.main %in% c("DC", "Macrophage", "Monocyte", "Neutrophils")]
Mye.SingleR.anno <- list(
  "main" = as.data.frame(SingleR(test = Mye.seu@assays$RNA@data, 
                                 ref = ref, 
                                 labels = ref$label.main)),
  "fine" = as.data.frame(SingleR(test = Mye.seu@assays$RNA@data, 
                                 ref = ref, 
                                 labels = ref$label.fine))
)

Mye.seu$SingleR.main.anno <- Mye.SingleR.anno$main$pruned.labels
Mye.seu$SingleR.fine.anno <- Mye.SingleR.anno$fine$pruned.labels
DimPlot(Mye.seu, group.by = "SingleR.main.anno", reduction = "umap_harmony")
DimPlot(Mye.seu, group.by = "SingleR.fine.anno", reduction = "umap_harmony")
saveRDS(Mye.seu, file.path(res.path, "Mye.seu.rds"))

# Switch label ------------------------------------------------------------

label <- data.frame("from" = unique(seu$Cell_type1))
label$to1 <- rep("unknown", nrow(label))
label$to2 <- rep("unknown", nrow(label))
write.table(label, file.path(res.path, "label.txt"),
            row.names = F, col.names = T, quote = F, sep = "\t")


# Fetch Annotation --------------------------------------------------------

LM22.annot <- function(gene){
  tmp <- unlist(lapply(LM22.signature, function(x) match(gene, x)))
  tmp <- tmp[!is.na(tmp)]
  if (length(tmp) == 0){
    return("")
  }else return(names(sort(tmp)[1]))
}

# Garbage -----------------------------------------------------------------

test <- RunFastMNN(object.list = SplitObject(seu, split.by = "Channel10X"))
features <- SelectIntegrationFeatures(object.list = SplitObject(seu, split.by = "Channel10X"))
anchors <- FindIntegrationAnchors(object.list = SplitObject(seu, split.by = "Channel10X"), 
                                  anchor.features = features)
tmp <- IntegrateData(anchorset = anchors)
tmp <- ScaleData(tmp, verbose = FALSE)
tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
tmp <- RunUMAP(tmp, reduction = "pca", dims = 1:15)
tmp <- FindNeighbors(tmp, reduction = "pca", dims = 1:15)
tmp <- FindClusters(tmp, resolution = 0.5)
FeaturePlot(tmp, features = plot.marker)
test <- RunFastMNN(object.list = SplitObject(seu, split.by = "Channel10X"))
