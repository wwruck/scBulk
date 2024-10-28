# scBulk
code and data for deconvolution of kidney and brain organoid trancriptome data

## Seurat analysis of sc-RNAseq kidney dataset GSE202109

Load libraries

``` r
library(Seurat)
library(gplots)
library(scMayoMap)
library(ggplot2)
library(svglite)
library(preprocessCore)
library(sva)
library(InstaPrism)
library(plotrix)
```

Get matrix, barcode and features .gz files of all samples (except the
CD45) from NCBI GEO **GSE202109** , matrix files:
GSM6094655_HKB14_matrix.mtx.gz, GSM6094658_HKB17_matrix.mtx.gz,
GSM6094668_HKB32_matrix.mtx.gz GSM6094656_HKB15_matrix.mtx.gz,
GSM6094659_HKB18_matrix.mtx.gz, GSM6094669_HKB33_matrix.mtx.gz
GSM6094657_HKB16_matrix.mtx.gz, GSM6094667_HKB31_matrix.mtx.gz,
GSM6094670_HKB34_matrix.mtx.gz and corresponding barcodes and features
files. Go to the directory where you downloaded the files.

``` r
setwd("D:/seurat/GSE202109") # replace with your working directory
expression_matrix <- ReadMtx( # M35
  mtx = "GSM6094655_HKB14_matrix.mtx.gz", features = "GSM6094655_HKB14_features.tsv.gz",
  cells = "GSM6094655_HKB14_barcodes.tsv.gz"
)
seurat_R1 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R1) <- "RNA"
seurat_R1[["percent.mt"]] <- PercentageFeatureSet(seurat_R1, pattern = "^MT-")
seurat_R1 <- subset(seurat_R1, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx( # M65
  mtx = "GSM6094656_HKB15_matrix.mtx.gz", features = "GSM6094656_HKB15_features.tsv.gz",
  cells = "GSM6094656_HKB15_barcodes.tsv.gz"
)
seurat_R2 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R2) <- "RNA"
seurat_R2[["percent.mt"]] <- PercentageFeatureSet(seurat_R2, pattern = "^MT-")
seurat_R2 <- subset(seurat_R2, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx(
  mtx = "GSM6094657_HKB16_matrix.mtx.gz", features = "GSM6094657_HKB16_features.tsv.gz",
  cells = "GSM6094657_HKB16_barcodes.tsv.gz"
)
seurat_R3 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R3) <- "RNA"
seurat_R3[["percent.mt"]] <- PercentageFeatureSet(seurat_R3, pattern = "^MT-")
seurat_R3 <- subset(seurat_R3, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx(
  mtx = "GSM6094658_HKB17_matrix.mtx.gz", features = "GSM6094658_HKB17_features.tsv.gz",
  cells = "GSM6094658_HKB17_barcodes.tsv.gz"
)
seurat_R4 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R4) <- "RNA"
seurat_R4[["percent.mt"]] <- PercentageFeatureSet(seurat_R4, pattern = "^MT-")
seurat_R4 <- subset(seurat_R4, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx( # F64
  mtx = "GSM6094659_HKB18_matrix.mtx.gz", features = "GSM6094659_HKB18_features.tsv.gz",
  cells = "GSM6094659_HKB18_barcodes.tsv.gz"
)
seurat_R5 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R5) <- "RNA"
seurat_R5[["percent.mt"]] <- PercentageFeatureSet(seurat_R5, pattern = "^MT-")
seurat_R5 <- subset(seurat_R5, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx( # F50
  mtx = "GSM6094667_HKB31_matrix.mtx.gz", features = "GSM6094667_HKB31_features.tsv.gz",
  cells = "GSM6094667_HKB31_barcodes.tsv.gz"
)
seurat_R6 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R6) <- "RNA"
seurat_R6[["percent.mt"]] <- PercentageFeatureSet(seurat_R6, pattern = "^MT-")
seurat_R6 <- subset(seurat_R6, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx(
  mtx = "GSM6094668_HKB32_matrix.mtx.gz", features = "GSM6094668_HKB32_features.tsv.gz",
  cells = "GSM6094668_HKB32_barcodes.tsv.gz"
)
seurat_R7 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R7) <- "RNA"
seurat_R7[["percent.mt"]] <- PercentageFeatureSet(seurat_R7, pattern = "^MT-")
seurat_R7 <- subset(seurat_R7, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx(
  mtx = "GSM6094669_HKB33_matrix.mtx.gz", features = "GSM6094669_HKB33_features.tsv.gz",
  cells = "GSM6094669_HKB33_barcodes.tsv.gz"
)
seurat_R8 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R8) <- "RNA"
seurat_R8[["percent.mt"]] <- PercentageFeatureSet(seurat_R8, pattern = "^MT-")
seurat_R8 <- subset(seurat_R8, subset = nFeature_RNA > 200 & percent.mt < 50)

expression_matrix <- ReadMtx(
  mtx = "GSM6094670_HKB34_matrix.mtx.gz", features = "GSM6094670_HKB34_features.tsv.gz",
  cells = "GSM6094670_HKB34_barcodes.tsv.gz"
)
seurat_R9 <- CreateSeuratObject(counts = expression_matrix)
DefaultAssay(seurat_R9) <- "RNA"
seurat_R9[["percent.mt"]] <- PercentageFeatureSet(seurat_R9, pattern = "^MT-")
seurat_R9 <- subset(seurat_R9, subset = nFeature_RNA > 200 & percent.mt < 50)
rm(expression_matrix)
gc()
```

    ##             used  (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells   7151489 382.0   13554787  724.0  13554787  724.0
    ## Vcells 101492249 774.4  274056140 2090.9 274041335 2090.8

``` r
slist=list(seurat_R1,seurat_R2,seurat_R3,seurat_R4,seurat_R5,seurat_R6,seurat_R7,seurat_R8,seurat_R9)
```

Seurat analysis: normalize and identify variable features for each
dataset independently

``` r
slist <- lapply(X = slist, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x,verbose = FALSE)
    x <- RunPCA(x, npcs = 30, verbose = FALSE)
})
```

select features for integration and integrate with faster version with
rpca

``` r
features <- SelectIntegrationFeatures(object.list = slist)
s.anchors <- FindIntegrationAnchors(object.list = slist,reduction="rpca", anchor.features = features)
s.combined <- IntegrateData(anchorset = s.anchors)
s.combined[["RNA"]] <- JoinLayers(s.combined[["RNA"]])
```

Run the standard workflow for visualization and clustering

``` r
s.combined <- ScaleData(s.combined, verbose = FALSE)
s.combined <- RunPCA(s.combined, npcs = 30, verbose = FALSE)
s.combined <- RunUMAP(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindNeighbors(s.combined, reduction = "pca", dims = 1:30)
s.combined <- FindClusters(s.combined, resolution = 0.8)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 27136
    ## Number of edges: 1165953
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8530
    ## Number of communities: 26
    ## Elapsed time: 5 seconds

``` r
rm(seurat_R1,seurat_R2,seurat_R3,seurat_R4,seurat_R5,seurat_R6,seurat_R7,seurat_R8,seurat_R9)
gc()
```

    ##             used   (Mb) gc trigger  (Mb)   max used    (Mb)
    ## Ncells   7332279  391.6   13554787   724   13554787   724.0
    ## Vcells 656403242 5008.0 1699214613 12964 1698929177 12961.9

Find markers for all clusters

``` r
DefaultAssay(s.combined) <- "RNA"
all.markers <- FindAllMarkers(s.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file="markers_GSE202109.txt",sep="\t",quote=F,row.names=F)
```

You can save the image here with save.image(“GSE202109.RData”) what will
take quite a while and create a huge file but speeds up if you want to
more frequently work with the cluster analysis results. You can then
load the data via the command load(“GSE202109.RData”). Otherwise you can
just continue without saving and re-loading the image.
