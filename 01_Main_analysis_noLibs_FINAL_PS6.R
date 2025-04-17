# Upstream analysis of scRNA of human NK cells ex vivo: QC, normalization,
# scaling, clustering and dim-reduction

rm(list=ls())

library(dplyr)
library(parallel)
library(Seurat)
library(scater)
library(scran)
library(purrr)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(data.table)
library(stringr)

path <- "~/Downloads/20440_NK/Clonal_NK/data/raw_ungunzipped/"

# Enable parallelization via mclapply
ncores <- detectCores()-1

# Sample names
my_samples <- c("CMVpos2", "CMVpos3", "CMVpos4", "CMVneg1", "CMVneg2")

# Set color palettes
grey_scale <- brewer.pal(n= 9 ,name="Greys")
red_scale <- brewer.pal(n= 9,name="Reds")
colorscale <- c(grey_scale[3], red_scale[2:9])
clusterpal <- brewer.pal(8, "Dark2")
clusterpal <- c(clusterpal[1:3],clusterpal[5:8])
NKG2Cpal <- c("#40BAD5","#120136")
divpal <- brewer.pal(n = 9, name = "RdBu")
adaptivepal <- c("#880E4F",  "#C2185B",  "#E91E63",  "#F06292", "#F8BBD0", "#4A148C", "#7B1FA2", "#9C27B0", "#BA68C8")

viriscale <- viridis(9)

#colorscale <- plasma(9)


##### Load the datasets #####
CMVpos2.data <- Read10X(data.dir = paste0(path, "data_CMVpos2"), strip.suffix = T)
CMVpos3.data <- Read10X(data.dir = paste0(path, "data_CMVpos3"), strip.suffix = T)
CMVpos4.data <- Read10X(data.dir = paste0(path, "data_CMVpos4"), strip.suffix = T)
CMVneg1.data <- Read10X(data.dir = paste0(path, "data_CMVneg1"), strip.suffix = T)
CMVneg2.data <- Read10X(data.dir = paste0(path, "data_CMVneg2"), strip.suffix = T)

# Initialize the Seurat objects with the raw (non-normalized data).
# Keep all features expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected features
CMVpos2 <- CreateSeuratObject(counts = CMVpos2.data, min.cells = 3, min.features = 200, project = "CMVpos2")
CMVpos3 <- CreateSeuratObject(counts = CMVpos3.data, min.cells = 3, min.features = 200, project = "CMVpos3")
CMVpos4 <- CreateSeuratObject(counts = CMVpos4.data, min.cells = 3, min.features = 200, project = "CMVpos4")
CMVneg1 <- CreateSeuratObject(counts = CMVneg1.data, min.cells = 3, min.features = 200, project = "CMVneg1")
CMVneg2 <- CreateSeuratObject(counts = CMVneg2.data, min.cells = 3, min.features = 200, project = "CMVneg2")

rm(CMVpos2.data, CMVpos3.data, CMVpos4.data, CMVneg1.data, CMVneg2.data)
CMVpos2
CMVpos3
CMVpos4
CMVneg1
CMVneg2
#LibA
#LibB

##### QC #####
# store mitochondrial percentage in object meta data
CMVpos2 <- PercentageFeatureSet(CMVpos2, pattern = "^MT-", col.name = "percent.mt")
CMVpos3 <- PercentageFeatureSet(CMVpos3, pattern = "^MT-", col.name = "percent.mt")
CMVpos4 <- PercentageFeatureSet(CMVpos4, pattern = "^MT-", col.name = "percent.mt")
CMVneg1 <- PercentageFeatureSet(CMVneg1, pattern = "^MT-", col.name = "percent.mt")
CMVneg2 <- PercentageFeatureSet(CMVneg2, pattern = "^MT-", col.name = "percent.mt")
CMVneg2 <- PercentageFeatureSet(CMVneg2, pattern = "^MT-", col.name = "percent.mt")

# Inspect Metadata
VlnPlot(CMVpos2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(CMVpos2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVpos3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVpos4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVpos4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVneg1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVneg1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
CombinePlots(plots = list(plot1, plot2))

plot1 <- FeatureScatter(CMVneg2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CMVneg2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Filter based on QC metrics
CMVpos2 <- subset(CMVpos2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 6)
CMVpos3 <- subset(CMVpos3, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVpos4 <- subset(CMVpos4, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg1 <- subset(CMVneg1, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg2 <- subset(CMVneg2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)

# Put Seurat objects into a list
my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)

# Visualize markers of contaminating cells
# my_SeuratList %>%
#   map(function(x) VlnPlot(object = x, features = c("HBA1", "HBA2", "IGJ"), slot = "counts"))

# # Filter erythrocytes and b cells out
# my_SeuratList %>%
#   map(function(x) subset(x, subset= HBA1 == 0 & HBA2 == 0 & IGJ == 0)) -> my_SeuratList
# my_SeuratList

my_SeuratList %>%
  map(~ VlnPlot(.x, features = c("HBA1", "HBA2", "IGJ"), layer = "counts"))

my_SeuratList <- map(my_SeuratList, function(x) {
  expr <- FetchData(x, vars = c("HBA1", "HBA2", "IGJ"), layer = "counts")
  keep_cells <- rowSums(expr) == 0
  subset(x, cells = rownames(expr)[keep_cells])
})
  
scran_normalize <- function(seu) {
    # Convert Seurat object to SingleCellExperiment
    sce <- as.SingleCellExperiment(seu)
    
    # Perform clustering for the sum factor calculation
    clusters <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters = clusters)
    
    # Print summary of size factors
    summary(sizeFactors(sce))
    
    # Normalize counts (log-transformed)
    sce <- logNormCounts(sce, log = FALSE)
    
    # Get the normalized counts matrix (normcounts)
    lognormcounts <- sce@assays@data$normcounts
    lognormcounts <- as(lognormcounts, "dgCMatrix")  # Convert to sparse matrix
    
    # Ensure dimnames match Seurat object
    dimnames(lognormcounts) <- dimnames(seu[["RNA"]]$counts)
    
    # Extract size factors and add to Seurat object metadata
    sizefactors <- sizeFactors(sce)
    names(sizefactors) <- names(seu$nCount_RNA)
    seu <- AddMetaData(seu, sizefactors, col.name = "sizefactors")
    
    # Store normalized counts in the RNA assay under the "counts" slot (not "data")
    seu[["RNA"]]$counts <- lognormcounts
    
    # Clean up
    rm(sce)
    
    return(seu)
  }

mclapply(my_SeuratList, mc.cores = ncores, function(x) scran_normalize(x)) -> my_SeuratList

# Visualize sizefactors
FeatureScatter(my_SeuratList[[2]], feature1 = "sizefactors", feature2 = "nCount_RNA")
FeatureScatter(my_SeuratList[[2]], feature1 = "sizefactors", feature2 = "nFeature_RNA")

##### Variable Features, Scaling #####
# Find variable features
my_SeuratList %>% 
  lapply(function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)) -> my_SeuratList

# Scale Data
lapply(my_SeuratList, function(x) NormalizeData(x)) -> my_SeuratList
lapply(my_SeuratList, function(x) ScaleData(x, features = rownames(x))) -> my_SeuratList

# Run PCA
lapply(my_SeuratList, function(x) RunPCA(x, verbose = FALSE)) -> my_SeuratList

# Inspect PCs in Elbowplot and DimHeatmaps
pdf(file = "PCA_results.pdf", title = "PCA Results", paper = "a4")
my_SeuratList %>% 
  map(function(x){
    DimHeatmap(x, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 19:27, cells = 500, balanced = TRUE, nfeatures = 22)
  })
my_SeuratList %>% 
  map(function(x) ElbowPlot(x))
dev.off()

# RunUMAP
my_SeuratList[[1]] <- RunUMAP(my_SeuratList[[1]], dims = 1:13, verbose = FALSE)
my_SeuratList[[2]] <- RunUMAP(my_SeuratList[[2]], dims = 1:14, verbose = FALSE)
my_SeuratList[[3]] <- RunUMAP(my_SeuratList[[3]], dims = 1:13, verbose = FALSE)
my_SeuratList[[4]] <- RunUMAP(my_SeuratList[[4]], dims = 1:10, verbose = FALSE)
my_SeuratList[[5]] <- RunUMAP(my_SeuratList[[5]], dims = 1:10, verbose = FALSE)
#my_SeuratList[[6]] <- RunUMAP(my_SeuratList[[6]], dims = 1:16, verbose = FALSE)
#my_SeuratList[[7]] <- RunUMAP(my_SeuratList[[7]], dims = 1:12, verbose = FALSE)

##### Clustering #####
# High resolution in some donors (i.e. overclustering)
# to identify contaminating ILCs (are in proximity to bright)
my_SeuratList[[1]] <- FindNeighbors(my_SeuratList[[1]], dims = 1:13, verbose = FALSE)
my_SeuratList[[1]] <- FindClusters(my_SeuratList[[1]], verbose = FALSE, resolution = 2, n.start = 100)

my_SeuratList[[2]] <- FindNeighbors(my_SeuratList[[2]], dims = 1:14, verbose = FALSE)
my_SeuratList[[2]] <- FindClusters(my_SeuratList[[2]], verbose = FALSE, resolution = 1, n.start = 100)

my_SeuratList[[3]] <- FindNeighbors(my_SeuratList[[3]], dims = 1:13, verbose = FALSE)
my_SeuratList[[3]] <- FindClusters(my_SeuratList[[3]], verbose = FALSE, resolution = 0.6, n.start = 100)

my_SeuratList[[4]] <- FindNeighbors(my_SeuratList[[4]], dims = 1:10, verbose = FALSE)
my_SeuratList[[4]] <- FindClusters(my_SeuratList[[4]], verbose = FALSE, resolution = 0.6, n.start = 100)

my_SeuratList[[5]] <- FindNeighbors(my_SeuratList[[5]], dims = 1:10, verbose = FALSE)
my_SeuratList[[5]] <- FindClusters(my_SeuratList[[5]], verbose = FALSE, resolution = 0.4, n.start = 100)

# Plot UMAPs with clustering and expression of ILC/NK markers
pdf(file = "UMAP_all_donors_ILCs.pdf", paper = "a4")
my_SeuratList %>% 
  map(function(x) DimPlot(x, label = T))
my_SeuratList %>% 
  map(function(x) FeaturePlot(x, features = c("CD56-adt", "IL7R", "GATA3", "IL2RA", "CD40LG"), cols = colorscale,
                              min.cutoff = "q1", max.cutoff = "q99", order = TRUE))
dev.off()



##### Removal of contaminating cells (mostly ILC2s: High IL7R, GATA3, IL2RA, CD40LG, negative for CD56) #####
my_SeuratList[[1]] <- subset(my_SeuratList[[1]], idents = "16", invert = TRUE) # Cluster 16 are ILC2
#my_SeuratList[[2]] <- subset(my_SeuratList[[2]], idents = "9", invert = TRUE) # Cluster 9 are ILC2
#my_SeuratList[[3]] <- subset(my_SeuratList[[3]], idents = "6", invert = TRUE) # Cluster 6 are ILC2
#my_SeuratList[[4]] <- subset(my_SeuratList[[4]], idents = "7", invert = TRUE) # Cluster 7 are ILC2
#my_SeuratList[[5]] <- subset(my_SeuratList[[5]], idents = "6", invert = TRUE) # Cluster 6 are ILC2
#my_SeuratList[[6]] <- subset(my_SeuratList[[6]], idents = "7", invert = TRUE) # Cluster 7 are ILC2
#my_SeuratList[[7]] <- subset(my_SeuratList[[7]], idents = "8", invert = TRUE) # Cluster 8 are ILC2


##### Re-run Dimred and Clustering #####
# Run PCA
my_SeuratList %>% 
  lapply(function(x) RunPCA(x, verbose = FALSE)) -> my_SeuratList

# Inspect PCs in Elbowplot and DimHeatmaps
pdf(file = "PCA_results_no_ILCs.pdf", title = "PCA Results", paper = "a4")
my_SeuratList %>% 
  map(function(x){
    DimHeatmap(x, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 22)
    DimHeatmap(x, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 22)
  })
my_SeuratList %>% 
  map(function(x) ElbowPlot(x))
dev.off()

# RunUMAP
my_SeuratList[[1]] <- RunUMAP(my_SeuratList[[1]], dims = 1:11, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[2]] <- RunUMAP(my_SeuratList[[2]], dims = 1:14, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[3]] <- RunUMAP(my_SeuratList[[3]], dims = 1:15, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[4]] <- RunUMAP(my_SeuratList[[4]], dims = 1:10, verbose = FALSE, n.neighbors = 20)
my_SeuratList[[5]] <- RunUMAP(my_SeuratList[[5]], dims = 1:13, verbose = FALSE, n.neighbors = 20)
# my_SeuratList[[6]] <- RunUMAP(my_SeuratList[[6]], dims = 1:14, verbose = FALSE, n.neighbors = 20)
# my_SeuratList[[7]] <- RunUMAP(my_SeuratList[[7]], dims = 1:14, verbose = FALSE, n.neighbors = 20)

#remove 6 and 7 (lib A and B) because we dont want those dang patients :)
my_SeuratList <- my_SeuratList[-c(6, 7)]

# Add donor annotation
CMVpos2$donor <- "CMVpos2"
CMVpos3$donor <- "CMVpos3"
CMVpos4$donor <- "CMVpos4"
CMVneg1$donor <- "CMVneg1"
CMVneg2$donor <- "CMVneg2"
#CMVpos1_5$donor <- "CMVpos1"
#CMVpos5_v3$donor <- "CMVpos5"


# Put Seurat objects into a list
my_SeuratList <- list(CMVpos2, CMVpos3, CMVneg1)
# my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)
# my_SeuratList <- list(CMVpos2, CMVpos3)
#rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)

# Assign serostatus
my_SeuratList[[1]]$serostatus <- "CMVpos"
my_SeuratList[[2]]$serostatus <- "CMVpos"
my_SeuratList[[3]]$serostatus <- "CMVneg" # hard coded this line FYI !!!!
# my_SeuratList[[4]]$serostatus <- "CMVneg"
# my_SeuratList[[5]]$serostatus <- "CMVneg"
# my_SeuratList[[6]]$serostatus <- "CMVpos"
# my_SeuratList[[7]]$serostatus <- "CMVpos"


##### Integration of datasets #####
# Find variable features, scale
library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize = 30000 * 1024^2) # for 128 Gb RAM

# Find overlapping variable features of CMVpos
shared_variable_features_CMVpos <- Reduce(intersect, list(VariableFeatures(my_SeuratList[[1]]),
                       VariableFeatures(my_SeuratList[[2]])#,
                       #VariableFeatures(my_SeuratList[[3]])
                       ))

shared_features <- Reduce(intersect, list(rownames(my_SeuratList[[1]]),
                                          rownames(my_SeuratList[[2]]),
                                          rownames(my_SeuratList[[3]])#,
                                          #rownames(my_SeuratList[[4]]),
                                          #rownames(my_SeuratList[[5]])
                                          ))

all_features <- Reduce(unique, list(rownames(my_SeuratList[[1]]),
                                       rownames(my_SeuratList[[2]]),
                                       rownames(my_SeuratList[[3]])#,
                                       #rownames(my_SeuratList[[4]]),
                                       #rownames(my_SeuratList[[5]])
                                    ))

# normalize counts and create data layer
my_SeuratList <- lapply(my_SeuratList, function(obj) {
  obj <- SketchData(
    object = obj,
    ncells = 2000,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  obj <- NormalizeData(obj, assay = "RNA")  # Apply normalization to the RNA assay
  return(obj)
})

my_SeuratList <- lapply(my_SeuratList, function(obj) {
  # Compute variable features if not already done
  if(length(VariableFeatures(obj)) == 0) {
    obj <- FindVariableFeatures(obj, assay = "RNA")  # Assure you're working with the RNA assay
  }
  return(obj)
})

# Integrate CMVpos
integration_anchors_CMVpos <- FindIntegrationAnchors(c(my_SeuratList[[1]], my_SeuratList[[2]]),
                                                     dims = 1:10, anchor.features = shared_variable_features_CMVpos, assay = c("RNA", "RNA"))
CMVpos <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
DefaultAssay(CMVpos) <- "integrated"

# integration_anchors_CMVpos <- FindIntegrationAnchors(c(my_SeuratList[[1]], my_SeuratList[[2]]),
#                                                      dims = 1:10, anchor.features = shared_variable_features_CMVpos, assay = c("RNA", "RNA"))
# CMVpos <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
# DefaultAssay(CMVpos) <- "integrated"

# Process integrated dataset
CMVpos <- ScaleData(CMVpos, vars.to.regress = "nCount_RNA")

CMVpos <- FindVariableFeatures(CMVpos, assay = "integrated") 

CMVpos <- RunPCA(CMVpos)

CMVpos <- RunUMAP(CMVpos, dims = 1:10, n.neighbors = 20, min.dist = 0.2)
CMVpos <- FindNeighbors(CMVpos, dims = 1:10)
CMVpos <- FindClusters(CMVpos, resolution = 0.2)

# Subcluster to resolve early dim
# temp <- subset(CMVpos, idents = "2")
# temp <- FindNeighbors(temp, dims = 1:10)
# temp <- FindClusters(temp, resolution = 0.2)
# DimPlot(temp)
# Idents(CMVpos, cells = WhichCells(temp, idents = "0")) <- "EarlyCD56dim"

DimPlot(CMVpos)
DimPlot(CMVpos, group.by = "orig.ident")

# chromatin accesibility (gene score)
# gene expression
FeaturePlot(CMVpos, features = c("KLRC2","ZBTB16", "GZMK","MKI67", "STMN1"))
Features(CMVpos)

CMVpos <- RenameIdents(CMVpos, "2" = "CD56bright", "1" = "CD56dim", "0" = "Adaptive", "3" = "Proliferating")
#levels(CMVpos) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating",  "Adaptive")
CMVpos$annotation <- Idents(CMVpos)

DimPlot(CMVpos, cols = c(clusterpal[1:4], adaptivepal))&NoLegend()

DimPlot(CMVpos, cols = c(clusterpal[1:4], adaptivepal), group.by = "ident") + 
  ggtitle("NK Cell Cluster Identity")


#export data

# for integrated pos data set
Idents(CMVpos) <- "annotation"
table(Idents(CMVpos))

markers_adaptive_bright <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "CD56bright")
write.csv(markers, file='CMVpos_adaptive_vs_bright.csv')

markers_adaptive_dim <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "CD56dim")
write.csv(markers, file='CMVpos_adaptive_vs_dim.csv')

markers_dim_adaptive <- FindMarkers(CMVpos, ident.1 = "CD56dim", ident.2 = "Adaptive")
write.csv(markers, file='CMVpos_dim_vs_adaptive.csv')

markers_adaptive_proliferating <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "Proliferating")
write.csv(markers, file='CMVpos_adaptive_vs_proliferating.csv')

markers_adaptive_bright <- FindMarkers(CMVpos, ident.1 = "CD56bright", ident.2 = "Adaptive")
write.csv(markers, file='CMVpos_bright_vs_adaptive.csv')

# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- markers_bright_adaptive[!is.na(markers_bright_adaptive$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)

# Install if needed
BiocManager::install("msigdbr")
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

BiocManager::install("fgsea")
library(fgsea)

BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
                      )

head(gsea_results[order(padj), ])


####

# for neg
head(CMVneg)
CMVneg@meta.data
Idents(CMVneg) <- "annotation"
markers <- FindMarkers(CMVneg, ident.1 = "CD56dim", ident.2 = "CD56bright")

# Process integrated dataset
CMVneg <- CMVneg1
DefaultAssay(CMVneg) <- "RNA"
CMVneg <- NormalizeData(CMVneg)
CMVneg <- ScaleData(CMVneg, vars.to.regress = "nCount_RNA")

CMVneg <- FindVariableFeatures(CMVneg, assay = "RNA") 
CMVneg <- RunPCA(CMVneg)
CMVneg <- RunUMAP(CMVneg, dims = 1:10)
DimPlot(CMVneg)

CMVneg <- FindNeighbors(CMVneg, dims = 1:10)
CMVneg <- FindClusters(CMVneg, resolution = 0.05, k.param=4)
DimPlot(CMVneg)

# gene expression
FeaturePlot(CMVneg, features = c("CX3CR1","GZMK","STMN1"))
Features(CMVneg)

# # Subcluster to resolve early dim
# temp <- subset(CMVneg, idents = "1")
# temp <- FindNeighbors(temp, dims = 1:10)
# temp <- FindClusters(temp, resolution = 0.2)
# DimPlot(temp)
# Idents(CMVneg, cells = WhichCells(temp, idents = "1")) <- "EarlyCD56dim"
# DimPlot(CMVneg)

CMVneg <- RenameIdents(CMVneg, "1" = "CD56bright", "0" = "CD56dim", "2" = "Proliferating")
#levels(CMVneg) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating")
CMVneg$annotation <- Idents(CMVneg)
DimPlot(CMVneg, group.by = "ident") +
  ggtitle("CMVneg NK Cell Cluster Identity")

# Assign populations also here
CMVneg$population <- "NKG2Cneg"
CMVneg$population[CMVneg$HTO_maxID=="NKG2Cpos"] <- "NKG2Cpos"


DimPlot(CMVneg, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE)

# Normalize RNA assay for full dataset
DefaultAssay(CMVpos) <- "RNA"
DefaultAssay(CMVneg) <- "RNA"
CMVpos <- NormalizeData(CMVpos)
CMVneg <- NormalizeData(CMVneg)
CMVpos <- ScaleData(CMVpos)
CMVneg <- ScaleData(CMVneg)

# Export integrated datasets
saveRDS(CMVpos, file = "CMVpos_scRNA")
saveRDS(CMVneg, file = "CMVneg_scRNA")

# Merge the integrated datasets
integrated_scRNA <- merge(CMVpos, y = CMVneg)
DefaultAssay(integrated_scRNA) <- "integrated"

# Scale, regress out depth effects
VariableFeatures(integrated_scRNA) <- intersect(VariableFeatures(CMVpos), rownames(CMVneg))

integrated_scRNA <- ScaleData(integrated_scRNA, features = VariableFeatures(integrated_scRNA), vars.to.regress = "nCount_RNA")

integrated_scRNA <- RunPCA(integrated_scRNA)

integrated_scRNA <- RunUMAP(integrated_scRNA, dims = c(1:10), min.dist = 0.2)

integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = c(1:10))
integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 0.2)
DimPlot(integrated_scRNA, cols = c(clusterpal[1:4], adaptivepal), pt.size = 0.1)

# Subcluster to resolve early CD56dim
temp <- subset(integrated_scRNA, idents = "2")
temp <- FindNeighbors(temp, dims = c(1:10))
temp <- FindClusters(temp, resolution = 0.2)
DimPlot(temp, pt.size = 0.8, label = T)
Idents(integrated_scRNA, cells =  WhichCells(temp, idents = "1")) <- "EarlyCD56dim"

integrated_scRNA <- RenameIdents(integrated_scRNA, "2" = "CD56bright", "0" = "CD56dim", "1" = "Adaptive", "3" = "Proliferating")
levels(integrated_scRNA) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating",  "Adaptive")
integrated_scRNA$annotation <- Idents(integrated_scRNA)

DimPlot(integrated_scRNA, cols = c(clusterpal[1:4], adaptivepal))


DimPlot(integrated_scRNA, group.by = "serostatus", cols = c("#47BCFF","#FF8A47"), pt.size = 0.6, shuffle = T)
DimPlot(integrated_scRNA, group.by = "serostatus", cols = c("#47BCFF","#FF8A47"), pt.size = 0.6, shuffle = T)+NoLegend()


# For downstream analysis, switch back to RNA assay
# and re-normalize to avoid depth effects
DefaultAssay(integrated_scRNA) <- "RNA"
integrated_scRNA <- NormalizeData(integrated_scRNA)
integrated_scRNA <- ScaleData(integrated_scRNA)

# Export processed dataset
saveRDS(integrated_scRNA, file = "integrated_scRNA")






