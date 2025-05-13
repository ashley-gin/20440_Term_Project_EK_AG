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
#LibA.data <- Read10X(data.dir = paste0(path, "data_CMVpos1_5/LibA/filtered_feature_bc_matrix"), strip.suffix = T)
#LibB.data <- Read10X(data.dir = paste0(path, "data_CMVpos1_5/LibB/filtered_feature_bc_matrix"), strip.suffix = T)

# Add suffix '-2' to cellnames of LibB to allow merging
#colnames(LibB.data) <- paste0(colnames(LibB.data), "-2")

# Initialize the Seurat objects with the raw (non-normalized data).
# Keep all features expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected features
CMVpos2 <- CreateSeuratObject(counts = CMVpos2.data, min.cells = 3, min.features = 200, project = "CMVpos2")
CMVpos3 <- CreateSeuratObject(counts = CMVpos3.data, min.cells = 3, min.features = 200, project = "CMVpos3")
CMVpos4 <- CreateSeuratObject(counts = CMVpos4.data, min.cells = 3, min.features = 200, project = "CMVpos4")
CMVneg1 <- CreateSeuratObject(counts = CMVneg1.data, min.cells = 3, min.features = 200, project = "CMVneg1")
CMVneg2 <- CreateSeuratObject(counts = CMVneg2.data, min.cells = 3, min.features = 200, project = "CMVneg2")
#LibA <- CreateSeuratObject(counts = LibA.data, min.cells = 3, min.features = 200, project = "CMVpos1_5")
#LibB <- CreateSeuratObject(counts = LibB.data, min.cells = 3, min.features = 200, project = "CMVpos1_5")

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
#LibA <- PercentageFeatureSet(LibA, pattern = "^MT-", col.name = "percent.mt")
#LibB <- PercentageFeatureSet(LibB, pattern = "^MT-", col.name = "percent.mt")

# Inspect Metadata
VlnPlot(CMVpos2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVpos4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(CMVneg2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(LibA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(LibB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


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

# plot1 <- FeatureScatter(LibA, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(LibA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# 
# plot1 <- FeatureScatter(LibB, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(LibB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))


# Filter based on QC metrics
CMVpos2 <- subset(CMVpos2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 6)
CMVpos3 <- subset(CMVpos3, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVpos4 <- subset(CMVpos4, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg1 <- subset(CMVneg1, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
CMVneg2 <- subset(CMVneg2, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
# LibA <- subset(LibA, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)
# LibB <- subset(LibB, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA < 10000 & percent.mt < 10)


# Put Seurat objects into a list
my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)
#my_SeuratList <- list(CMVpos2, CMVpos3)
#my_SeuratList <- list(LibA, LibB)
#rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2)

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


#Unable to combine ATAC-seq data due to unannotated ADT labels
# Add CITE-Seq data for LibA, LibB separately,for those it was counted with kite, so take these out for now
# LibA <- my_SeuratList[[6]]
# LibB <- my_SeuratList[[7]]
# 
# CMVpos1_5 <- merge(LibA, LibB, merge.data = TRUE, project = "CMVpos1_5")
# 
# my_SeuratList <- c(my_SeuratList[[1]], my_SeuratList[[2]], my_SeuratList[[3]],
#                    my_SeuratList[[4]], my_SeuratList[[5]])


##### Add normalized CITE-Seq data and exclude doublets #####
# Define marker/hashtag panel and assign populations in the same order
# 
# my_markers <- c("CD62L-adt", "CD57-adt", "CD161-adt", "Anti_Biotin_NKG2A-adt", "CD56-adt",
#                 "CD16-adt", "CD2-adt", "KLRG1-adt", "CXCR4-adt", "KIR2DL2_L3-adt", "KIR2DL1_S135-adt",
#                 "KIR3DL1-adt", "LAG3-adt", "NKp30-adt")
# my_hashtags <- c("Hashtag1-adt", "Hashtag2-adt")
# my_populations <- c("NKG2Cpos", "NKG2Cneg") 
# 
# Add_normCITEseq <- function(my_Seurat, name = character()){
#   
#   # Read in ADT/HTO data
#   my_ADT <- as.sparse((x = read.delim(file = paste0(paste0(path, "data_"), name,
#                                                     "/CITE-Seq/ATCACGAT.umi.txt"),
#                                       header = TRUE, row.names = 1)))
#   my_ADT <- my_ADT[rownames(my_ADT)%in% my_markers,] 
#   my_HTO <- as.sparse((x = read.delim(file = paste0(paste0(path, "data_"), name,
#                                                     "/CITE-Seq/ATTACTCG.umi.txt"),
#                                       header = TRUE, row.names = 1)))
#   my_HTO <- my_HTO[rownames(my_HTO)%in% my_hashtags,] 
#   my_HTO
#   
#   # Define HTOs
#   rownames(my_HTO) <- my_populations
#   
#   # Trim data to cells present in transcriptomes and ADT/HTO
#   cellfilter.adt <- colnames(my_ADT) %in% colnames(my_Seurat)
#   cellfilter.transcriptomes <- colnames(my_Seurat)[colnames(my_Seurat) %in% colnames(my_ADT)]
#   my_ADT <- my_ADT[,cellfilter.adt]
#   my_Seurat <- subset(my_Seurat, cells = cellfilter.transcriptomes)
#   cellfilter.hto <- colnames(my_HTO) %in% colnames(my_Seurat)
#   my_HTO <- my_HTO[,cellfilter.hto]
#   
#   # Add ADT/HTO to Seurat
#   my_Seurat[["ADT"]] <- CreateAssayObject(counts = my_ADT)
#   my_Seurat[["HTO"]] <- CreateAssayObject(counts = my_HTO)
#   
#   # Normalize
#   my_Seurat <- NormalizeData(object = my_Seurat, assay = "ADT", normalization.method = "CLR")
#   my_Seurat <- NormalizeData(object = my_Seurat, assay = "HTO", normalization.method = "CLR")
#   my_Seurat <- ScaleData(object = my_Seurat, assay = "ADT")
#   
#   # Remove doublets
#   my_Seurat <- HTODemux(my_Seurat, assay = "HTO", positive.quantile = 0.99)
#   print(table(my_Seurat$HTO_classification.global))
#   my_Seurat <- my_Seurat[,my_Seurat$HTO_classification.global=="Singlet"]
#   
#   return(my_Seurat)
# } 
# 
# my_SeuratList %>% 
#   map2(my_samples[1:5], function(x, y) Add_normCITEseq(x, name = y)) -> my_SeuratList
# 
# ##### Addition of CITE-Seq data quantified with bustools to LibA/LibB#####
# # Import data
# import_kite_counts <- function(path){
#   mtx <- fread(paste0(path,"featurecounts.mtx"), header = FALSE)
#   dim <- mtx[1,]
#   mtx <- mtx[-1,]
#   matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
#   rownames(matx) <- fread(paste0(path,"featurecounts.barcodes.txt"), header = FALSE)[[1]]
#   colnames(matx) <- (fread(paste0(path,"featurecounts.genes.txt"), header = FALSE)[[1]])
#   return(t(matx))
# }
# 
# LibA_HTO <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibA_HTO/featurecounts/"))
# LibB_HTO <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibB_HTO/featurecounts/"))
# 
# LibA_ADT <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibA_ADT/featurecounts/"))
# LibB_ADT <- import_kite_counts(paste0(path, "data_CMVpos1_5/CMVpos1_5_CITE/LibB_ADT/featurecounts/"))
# 
# 
# 
# 
# # Harmonize cellnames
# colnames(LibA_HTO) <- paste0(colnames(LibA_HTO), "-1")
# table(colnames(LibA_HTO)%in%colnames(LibA))
# LibA_HTO <- LibA_HTO[,colnames(LibA_HTO)%in%colnames(LibA)]
# 
# colnames(LibA_ADT) <- paste0(colnames(LibA_ADT), "-1")
# table(colnames(LibA_ADT)%in%colnames(LibA))
# LibA_ADT <- LibA_ADT[,colnames(LibA_ADT)%in%colnames(LibA)]
# 
# colnames(LibB_HTO) <- paste0(colnames(LibB_HTO), "-2")
# table(colnames(LibB_HTO)%in%colnames(LibB))
# LibB_HTO <- LibB_HTO[,colnames(LibB_HTO)%in%colnames(LibB)]
# 
# colnames(LibB_ADT) <- paste0(colnames(LibB_ADT), "-2")
# table(colnames(LibB_ADT)%in%colnames(LibB))
# LibB_ADT <- LibB_ADT[,colnames(LibB_ADT)%in%colnames(LibB)]
# 
# # Harmonize rownames with other objects
# rownames(LibA_ADT) <- paste0(rownames(LibA_ADT), "-adt")
# rownames(LibA_ADT)[rownames(LibA_ADT) == "anti-Biotin-adt"] <- "Anti-Biotin-NKG2A-adt"
# rownames(LibA_ADT)[rownames(LibA_ADT) == "KIR2DL1-S1-S3-S5-adt"] <- "KIR2DL1-S135-adt"
# rownames(LibA_ADT)[rownames(LibA_ADT) == "LAG-3-adt"] <- "LAG3-adt"
# 
# rownames(LibB_ADT) <- paste0(rownames(LibB_ADT), "-adt")
# rownames(LibB_ADT)[rownames(LibB_ADT) == "anti-Biotin-adt"] <- "Anti-Biotin-NKG2A-adt"
# rownames(LibB_ADT)[rownames(LibB_ADT) == "KIR2DL1-S1-S3-S5-adt"] <- "KIR2DL1-S135-adt"
# rownames(LibB_ADT)[rownames(LibB_ADT) == "LAG-3-adt"] <- "LAG3-adt"
# 
# # One hashtag is mislabelled (CMVneg1 instead of CMVpos1), correct
# rownames(LibA_HTO)[4] <- "CMVpos1_NKG2Cneg"
# rownames(LibB_HTO)[4] <- "CMVpos1_NKG2Cneg"
# 
# 
# # Create Assay Objects
# LibA_HTO_assay <- CreateAssayObject(counts = LibA_HTO, min.features = 1)
# LibB_HTO_assay <- CreateAssayObject(counts = LibB_HTO, min.features = 1)
# 
# 
# 
# 
# # Merge
# CMVpos1_5_HTO <- CreateAssayObject(counts= cbind(LibA_HTO_assay@counts, LibB_HTO_assay@counts))
# 
# # Include only cells which are present in both modalities
# CMVpos1_5 <- subset(CMVpos1_5, cells = Cells(CMVpos1_5_HTO))
# 
# # Add HTO to objects
# CMVpos1_5[["HTO"]] <- CMVpos1_5_HTO
# 
# 
# 
# # There is a problem with hashtag 1 => Allmost all cells are marked by it and the higher counts
# # in the actual samples are not sufficient to get precise demultiplexing.
# # Luckily, all others work nicely, so we can correct for this by selecting all cells
# # which are negative for all other hashtags (thresholds set based on violin plots) and
# # setting the value for hashtag 1 to 100, all others to 0. This way we can still use HTODemux
# VlnPlot(CMVpos1_5, features = rownames(CMVpos1_5[["HTO"]][1:2]), log = TRUE, ncol =2)
# VlnPlot(CMVpos1_5, features = rownames(CMVpos1_5[["HTO"]][3:4]), log = TRUE, ncol =2)
# table(LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30)
# table(LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30)
# 
# # Add a normal distribution around 8 to be able to visualize
# LibA_HTO_corrected <- LibA_HTO
# LibB_HTO_corrected <- LibB_HTO
# 
# LibA_HTO_corrected[1,] <- round(rnorm(ncol(LibA_HTO), mean = 5, sd = 1))
# LibB_HTO_corrected[1,] <- round(rnorm(ncol(LibB_HTO), mean = 5, sd = 1))
# 
# LibA_HTO_corrected[1,LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30] <- round(rnorm(ncol(LibA_HTO_corrected[,LibA_HTO[2,]<10&LibA_HTO[3,]<30&LibA_HTO[4,]<30]), mean = 100, sd = 15))
# LibB_HTO_corrected[1,LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30] <- round(rnorm(ncol(LibB_HTO_corrected[,LibB_HTO[2,]<10&LibB_HTO[3,]<30&LibB_HTO[4,]<30]), mean = 100, sd = 15))
# 
# # Create assay
# # Create Assay Objects
# LibA_HTO_assay_corrected <- CreateAssayObject(counts = LibA_HTO_corrected, min.features = 1)
# LibB_HTO_assay_corrected <- CreateAssayObject(counts = LibB_HTO_corrected, min.features = 1)
# 
# # Merge 
# CMVpos1_5_HTO_corrected <- CreateAssayObject(counts= cbind(LibA_HTO_assay_corrected@counts, LibB_HTO_assay_corrected@counts))
# 
# # Add corrected assay to Seurat object
# CMVpos1_5[["HTO_corrected"]] <- CMVpos1_5_HTO_corrected
# 
# 
# 
# # Normalize
# CMVpos1_5 <- NormalizeData(CMVpos1_5, assay = "HTO_corrected", normalization.method = "CLR")
# 
# # Demultiplex
# CMVpos1_5<- HTODemux(CMVpos1_5, assay = "HTO_corrected", positive.quantile = 0.995,)
# table(CMVpos1_5$HTO_corrected_classification.global)
# CMVpos1_5 <- subset(CMVpos1_5, HTO_corrected_classification.global == "Singlet")
# 
# table(CMVpos1_5$HTO_corrected_maxID)
# 
# # Check how the distribution looks
# RidgePlot(CMVpos1_5, assay = "HTO_corrected", 
#           features = rownames(CMVpos1_5[["HTO_corrected"]]), ncol = 2, log = T)
# 
# 
# 
# # Add ADT data
# LibA_ADT <- CreateAssayObject(counts = LibA_ADT)
# LibB_ADT<- CreateAssayObject(counts = LibB_ADT)
# 
# 
# # Normalize per Library
# LibA_ADT <- NormalizeData(LibA_ADT, normalization.method = "CLR")
# LibB_ADT <- NormalizeData(LibB_ADT, normalization.method = "CLR")
# 
# 
# 
# 
# # Merge, scale and add to Objects
# CMVpos1_5_ADT <- CreateAssayObject(data=cbind(LibA_ADT@data, LibB_ADT@data))
# CMVpos1_5_ADT <- subset(CMVpos1_5_ADT, cells = Cells(CMVpos1_5))
# 
# 
# CMVpos1_5[["ADT"]] <- CMVpos1_5_ADT
# 
# # Demultiplex CMVpos1/5 from LibA/LibB
# CMVpos1_v3 <- subset(CMVpos1_5, `HTO_corrected_maxID`%in%c("CMVpos1-NKG2Cpos", "CMVpos1-NKG2Cneg"))
# CMVpos5_v3 <- subset(CMVpos1_5, `HTO_corrected_maxID`%in%c("CMVpos5-NKG2Cpos", "CMVpos5-NKG2Cneg"))
# 
# 
# my_SeuratList <- c(my_SeuratList, CMVpos1_v3, CMVpos5_v3)
# 
# # Export QCed datasets
# my_SeuratList %>% 
#   map2(my_samples, function(x, y) saveRDS(x, file = paste0(y, "_QC_ADT_HTO")))
# 
# # Import QCed datasets
# CMVpos2 <- readRDS(file = "CMVpos2_QC_ADT_HTO")
# CMVpos3 <- readRDS(file = "CMVpos3_QC_ADT_HTO")
# CMVpos4 <- readRDS(file = "CMVpos4_QC_ADT_HTO")
# CMVneg1 <- readRDS(file = "CMVneg1_QC_ADT_HTO")
# CMVneg2 <- readRDS(file = "CMVneg2_QC_ADT_HTO")
# CMVpos1_v3 <- readRDS(file = "CMVpos1_v3_QC_ADT_HTO")
# CMVpos5_v3 <- readRDS(file = "CMVpos5_v3_QC_ADT_HTO")
# 
# my_SeuratList <- list(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)
# rm(CMVpos2, CMVpos3, CMVpos4, CMVneg1, CMVneg2, CMVpos1_v3, CMVpos5_v3)
# 

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

#my_SeuratList[[6]] <- FindNeighbors(my_SeuratList[[6]], dims = 1:16, verbose = FALSE)
#my_SeuratList[[6]] <- FindClusters(my_SeuratList[[6]], verbose = FALSE, resolution = 0.6, n.start = 100)

#my_SeuratList[[7]] <- FindNeighbors(my_SeuratList[[7]], dims = 1:12, verbose = FALSE)
#my_SeuratList[[7]] <- FindClusters(my_SeuratList[[7]], verbose = FALSE, resolution = 0.6, n.start = 100)

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

#remove 6 and 7 (lib A and B)
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
CMVpos_int <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
DefaultAssay(CMVpos_int) <- "integrated"

# integration_anchors_CMVpos <- FindIntegrationAnchors(c(my_SeuratList[[1]], my_SeuratList[[2]]),
#                                                      dims = 1:10, anchor.features = shared_variable_features_CMVpos, assay = c("RNA", "RNA"))
# CMVpos <- IntegrateData(integration_anchors_CMVpos, dims = 1:10, features.to.integrate = all_features)
# DefaultAssay(CMVpos) <- "integrated"

# Process integrated dataset
CMVpos <- ScaleData(CMVpos_int, vars.to.regress = "nCount_RNA")

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

library(viridisLite) 
# chromatin accesibility (gene score)
# gene expression
FeaturePlot(CMVpos, features = c("KLRC2","ZBTB16", "GZMK","MKI67", "STMN1"),   cols = viridis_colors,
            min.cutoff=0,
            max.cutoff=4)
Features(CMVpos)

CMVpos <- RenameIdents(CMVpos, "2" = "CD56bright", "1" = "CD56dim", "0" = "Adaptive", "3" = "Proliferating")
#levels(CMVpos) <- c("CD56bright", "EarlyCD56dim","CD56dim", "Proliferating",  "Adaptive")
CMVpos$annotation <- Idents(CMVpos)

DimPlot(CMVpos, cols = c(clusterpal[1:4], adaptivepal))&NoLegend()

DimPlot(CMVpos, cols = c(clusterpal[1:4], adaptivepal), group.by = "ident") + 
  ggtitle("NK Cell Cluster Identity")

#####ADAPTIVE v. CD56Bright########

# export data comparing adaptive NK to other subtypes

# for integrated pos data set
Idents(CMVpos) <- "annotation"
table(Idents(CMVpos))

markers_adaptive_bright <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "CD56bright")
write.csv(markers_adaptive_bright, file='CMVpos_adaptive_vs_bright.csv')

markers_adaptive_dim <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "CD56dim")
write.csv(markers_adaptive_dim, file='CMVpos_adaptive_vs_dim.csv')

markers_adaptive_proliferating <- FindMarkers(CMVpos, ident.1 = "Adaptive", ident.2 = "Proliferating")
write.csv(markers_adaptive_proliferating, file='CMVpos_adaptive_vs_proliferating.csv')

# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- markers_adaptive_bright[!is.na(markers_adaptive_bright$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))
#ranks <- -ranks

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)

# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:10]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:10]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))

# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: Adaptive vs. CD56 Bright") +
      theme(
        axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 14, face = "bold")
  )

# plotting all (not just top 5 and bottom 5)
ggplot(gsea_plot_data, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       fill = "Significant", title = "Adaptive vs. CD56 Bright Pathway Enrichment") +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )


#########ADAPTIVE v. CD56dim###########

markers_clean <- markers_adaptive_dim[!is.na(markers_adaptive_dim$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))
#ranks <- -ranks

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)

# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:10]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:10]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))

# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: Adaptive vs. CD56 Dim"
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )



############# Adaptive vs. Proliferating ######


markers_clean <- markers_adaptive_proliferating[!is.na(markers_adaptive_proliferating$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))
#ranks <- -ranks

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)

# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:10]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:10]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))

# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: Adaptive vs. Proliferating"
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )




#######VOLCANO PLOT#############


### volcano plot !!!!!!

# set the Identity Class for Comparison

CMVpos_int$cmv_status <- "Positive"
CMVneg1$cmv_status <- "Negative"

# Merge into one Seurat object
combined_seurat <- merge(CMVpos_int, y = CMVneg1, add.cell.ids = c("CMVpos", "CMVneg1"))
Idents(combined_seurat) <- "cmv_status"

# Normalize the merged object

DefaultAssay(combined_seurat) <- "RNA"

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)

combined_seurat <- JoinLayers(combined_seurat)

markers_cmv <- FindMarkers(
  combined_seurat,
  ident.1 = "Positive",
  ident.2 = "Negative",
  logfc.threshold = 0,
  min.pct = 0.1,
  assay = "RNA"
)



markers_cmv$gene <- rownames(markers_cmv)
markers_cmv$logP <- -log10(markers_cmv$p_val_adj + 1e-300)
markers_cmv$significant <- with(markers_cmv, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

library(ggplot2)
library(dplyr)

top_up <- markers_cmv %>%
  filter(significant & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- markers_cmv %>%
  filter(significant & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

library(ggplot2)
library(ggrepel)  # For better text labels

library(ggplot2)
library(ggrepel)

markers_cmv$regulation <- "Not Significant"
markers_cmv$regulation[markers_cmv$avg_log2FC >= 1 ] <- "Upregulated"
markers_cmv$regulation[markers_cmv$avg_log2FC <= -1] <- "Downregulated"

ggplot(markers_cmv, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = regulation), alpha = 0.6, size=2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not Significant" = "gray")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 5,
                  fontface = "bold",
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in CMV+ vs CMV− Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")+
  theme(legend.position = "none")


## old
ggplot(markers_cmv, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in CMV+ vs CMV− Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +  # Adjust as needed
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")

install.packages("viridisLite")  # or viridis

# # ok now let's plot those 5 upregulated genes on the umap and see where they fall
# FeaturePlot(CMVpos, features = c("CD3D", "CD3E","MYOM2", "RPS26", "HLA-C", "KLRC2", "HLA-DPB1", "CCL5", "THEMIS2"), 
#             cols = c("red", "white", "blue")  # low → mid → high expression
# )

# # ok now let's plot those 5 upregulated genes on the umap and see where they fall
Features(CMVpos)
library(viridisLite)  # or library(viridis)
viridis_colors <- viridis(2)  # or more, e.g., viridis(10) for smoother gradient

FeaturePlot(
  CMVpos,
  features = c("CD3D", "CD3E", "MYOM2", "RPS26", "HLA-C", "KLRC2", "HLA-DPB1", "CCL5", "THEMIS2"),
  cols = viridis_colors,
  min.cutoff=0,
  max.cutoff=4
)

FeaturePlot(
  CMVpos,
  features = c("CD3D", "CD3E", "MYOM2", "RPS26", "HLA-C", "KLRC2", "HLA-DPB1", "CCL5", "THEMIS2")
)

# and downregulated genes - update genes here
FeaturePlot(CMVpos, 
            features = c("IFI44","IFI44L", "RP11-166O4.6","MX1","XCL1", "IFITM3", "ISG15", "LY6E", "FCER1G"),
            cols = viridis_colors,
            min.cutoff=0,
            max.cutoff=4
            )
Features(CMVpos)

# Install patchwork if you haven't already
install.packages("patchwork")
library(patchwork)

# Create individual plots but suppress legend on all except the first
plots <- lapply(c("CD3D", "CD3E", "MYOM2", "RPS26", "HLA-C", "KLRC2", "HLA-DPB1", "CCL5", "THEMIS2"), function(gene) {
  FeaturePlot(CMVpos, features = gene, cols = viridis(2)) + 
    theme(legend.position = "none")
})

# Add the legend from one selected gene (e.g. "IFI44")
legend_plot <- FeaturePlot(CMVpos, features = "CD3D", cols = viridis(2))

# Combine all with legend at the bottom
wrap_plots(plots) + plot_layout(guides = "collect") & theme(legend.position = "right")

# Create individual plots but suppress legend on all except the first
plots <- lapply(c("IFI44","IFI44L", "RP11-166O4.6","MX1","XCL1", "IFITM3", "ISG15", "LY6E", "FCER1G"), function(gene) {
  FeaturePlot(CMVpos, features = gene, cols = viridis(100)) + 
    theme(legend.position = "none")
})

# Add the legend from one selected gene (e.g. "IFI44")
legend_plot <- FeaturePlot(CMVpos, features = "IFI44", cols = viridis(100))

# Combine all with legend at the bottom
wrap_plots(plots) + plot_layout(guides = "collect") & theme(legend.position = "right")

# GSEA for memory-like phenotype

# Create a named vector of log2FC values, named by gene
ranked_genes <- markers_cmv %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC)) %>%  # or arrange(avg_log2FC) if needed
  tibble::deframe()

# Check the output
head(ranked_genes)

# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE) # elif 
ranked_genes <- -ranked_genes

# Check the output
head(ranked_genes)

library(fgsea)
#BiocManager::install("GSEABase")
library(GSEABase)

# Load the gene set
gmt_file <- "ImmuneSigDB.gmt"
geneset <- getGmt(gmt_file)

# Convert gene set
pathways <- geneIds(geneset)

#BiocManager::install("msigdbr")
library(msigdbr)

# Run fgsea
fgsea_res <- fgsea(pathways = pathways, stats = ranked_genes)

# View results
head(fgsea_res)

# Order results by NES
gsea_plot_goldrath <- as.data.table(fgsea_res)
gsea_plot_goldrath <- gsea_plot_goldrath[order(NES)]
gsea_plot_goldrath$gs_name <- factor(gsea_plot_goldrath$pathway, levels = gsea_plot_goldrath$pathway)

# create significant definition
gsea_plot_goldrath$significant <- factor(gsea_plot_goldrath$padj < 0.05, levels = c(TRUE, FALSE))

gsea_plot_goldrath_significant <- gsea_plot_goldrath[gsea_plot_goldrath$significant == TRUE, ]

filtered_res <- gsea_plot_goldrath_significant %>%
  filter(str_detect(tolower(pathway), "nk|pathway|cd8|memory"))

# Optional: sort by NES or p-value if you want better plot ordering
filtered_res <- filtered_res %>% arrange(desc(NES))


ggplot(filtered_res, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
  fill = "Significant", title = "Pathways upregulated in HCMV Positive vs Negative Patients") +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )

#GSEA of combined CMV+ and CMV- for Hallmark gene sets

# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- markers_cmv[!is.na(markers_cmv$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)
ranks <- -ranks
# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:10]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:10]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))


# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: All NK cells CMV+ v. CMV- Patients"
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )

upregulated <- rownames(markers_cmv[markers_cmv$avg_log2FC > 0, ])
top_genes <- intersect(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], upregulated)

top_genes <- intersect(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], rownames(markers_cmv))
DoHeatmap(cd56_dim, features = top_genes)

Idents(combined_seurat_subset) <- "cmv_status"
upregulated <- rownames(markers_cmv[markers_cmv$avg_log2FC > 0, ])
top_genes <- intersect(pathways[["GSE21360_NAIVE_VS_QUATERNARY_MEMORY_CD8_TCELL_DN"]], upregulated)
DoHeatmap(
  object = combined_seurat_subset,
  features = top_genes,
  assay = "integrated",             # Or whatever assay you're using
  slot = "data"              # Usually "data" for normalized expression
)



DoHeatmap(combined_seurat_subset, features = top_genes)

################

# Subset for only CD56dim cells and tell it to look at cmv_status ident to pull the positive vs negative info
cd56_dim <- subset(combined_seurat, subset = annotation == "CD56dim")
Idents(cd56_dim) <- "cmv_status"

# Differential expression
de_markers_cd56dim <- FindMarkers(cd56_dim, ident.1 = "Positive", ident.2 = "Negative", 
                                  logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
# Add gene names as a column
de_markers_cd56dim$gene <- rownames(de_markers)

# Classify significance and direction
de_markers$significance <- "Not significant"
de_markers$significance[de_markers$avg_log2FC > 0.25 & de_markers$p_val_adj < 0.05] <- "Upregulated"
de_markers$significance[de_markers$avg_log2FC < -0.25 & de_markers$p_val_adj < 0.05] <- "Downregulated"

# Factor levels to control legend order
de_markers$significance <- factor(de_markers$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

# Plot
ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", color = "Gene Regulation") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")



######### CD56 comparison between positive vs negative patients ######### 

# do a new combined seurat

# set the Identity Class for Comparison
CMVpos$cmv_status <- "Positive"
CMVneg$cmv_status <- "Negative"

# Merge into one Seurat object
combined_seurat_subtypes_1to1 <- merge(CMVpos, y = CMVneg, add.cell.ids = c("CMVpos", "CMVneg"))

DefaultAssay(combined_seurat_subtypes_1to1) <- "RNA"

combined_seurat_subtypes_1to1 <- NormalizeData(combined_seurat_subtypes_1to1)
combined_seurat_subtypes_1to1 <- FindVariableFeatures(combined_seurat_subtypes_1to1, nfeatures = 2000)

###
# Set a seed for reproducibility (optional)
set.seed(42)

# Randomly subset 5000 cells from your Seurat object
subset_cells <- sample(Cells(combined_seurat_subtypes_1to1), 5000)
combined_seurat_subset <- subset(combined_seurat_subtypes_1to1, cells = subset_cells)

combined_seurat_subset <- NormalizeData(combined_seurat_subset)
combined_seurat_subset <- FindVariableFeatures(combined_seurat_subset, nfeatures = 2000)
combined_seurat_subset <- JoinLayers(combined_seurat_subset)


######################  subtype CD56 dim ###################### 

cd56_dim <- subset(combined_seurat_subset, subset = annotation == "CD56dim")
Idents(cd56_dim) <- "cmv_status"

de_markers_cd56dim <- FindMarkers(cd56_dim, ident.1 = "Positive", ident.2 = "Negative", 
                                  logfc.threshold = 0, min.pct = 0.1, assay = "RNA")

# Add gene names as a column
de_markers_cd56dim$gene <- rownames(de_markers_cd56dim)

# Classify significance and direction
de_markers_cd56dim$significance <- "Not significant"
de_markers_cd56dim$significance[de_markers_cd56dim$avg_log2FC > 0.25 & de_markers_cd56dim$p_val_adj < 0.05] <- "Upregulated"
de_markers_cd56dim$significance[de_markers_cd56dim$avg_log2FC < -0.25 & de_markers_cd56dim$p_val_adj < 0.05] <- "Downregulated"

de_markers_cd56dim$gene <- rownames(de_markers_cd56dim)
de_markers_cd56dim$logP <- -log10(de_markers_cd56dim$p_val_adj + 1e-300)
de_markers_cd56dim$significant <- with(de_markers_cd56dim, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

de_markers_cd56dim$regulation <- "Not Significant"
de_markers_cd56dim$regulation[de_markers_cd56dim$avg_log2FC >= 1 ] <- "Upregulated"
de_markers_cd56dim$regulation[de_markers_cd56dim$avg_log2FC <= -1] <- "Downregulated"

# Factor levels to control legend order
de_markers_cd56dim$significance <- factor(de_markers_cd56dim$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

top_up <- de_markers_cd56dim %>%
  filter(significant & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- de_markers_cd56dim %>%
  filter(significant & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)


ggplot(de_markers_cd56dim, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = regulation), alpha = 0.6, size=2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not Significant" = "gray")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 5,
                  fontface = "bold",
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in CD56 Dim Subtype in CMV+ vs CMV- Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold")    # Axis tick labels
        )



###CD56 dim GSEA########
# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- de_markers_cd56dim[!is.na(de_markers_cd56dim$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)
ranks <- -ranks
# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:5]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:5]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))


# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: CD56 Dim in CMV+ v. CMV- Patients"
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Make y-axis labels bold and larger
    axis.text.x = element_text(size = 10),                    # Optional: adjust x-axis label size
    axis.title = element_text(size = 14, face = "bold"),      # Axis titles bold and larger
    plot.title = element_text(size = 16, face = "bold"),       # Title bold and larger
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )


upregulated <- rownames(de_markers_cd56dim[de_markers_cd56dim$avg_log2FC > 0, ])
upregulated <- -upregulated
top_genes <- intersect(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], upregulated)

top_genes <- intersect(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], rownames(de_markers_cd56dim))
DoHeatmap(cd56_dim, features = top_genes)


######################  subtype CD56 bright ###################### 

cd56_bright <- subset(combined_seurat_subset, subset = annotation == "CD56bright")
Idents(cd56_bright) <- "cmv_status"

de_markers_cd56bright <- FindMarkers(cd56_bright, ident.1 = "Positive", ident.2 = "Negative", 
                                  logfc.threshold = 0, min.pct = 0.1, assay = "RNA")

# Add gene names as a column
de_markers_cd56bright$gene <- rownames(de_markers_cd56bright)

# Classify significance and direction
de_markers_cd56bright$significance <- "Not significant"
de_markers_cd56bright$significance[de_markers_cd56bright$avg_log2FC > 0.25 & de_markers_cd56bright$p_val_adj < 0.05] <- "Upregulated"
de_markers_cd56bright$significance[de_markers_cd56bright$avg_log2FC < -0.25 & de_markers_cd56bright$p_val_adj < 0.05] <- "Downregulated"

# Factor levels to control legend order
de_markers_cd56bright$significance <- factor(de_markers_cd56bright$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

de_markers_cd56bright$gene <- rownames(de_markers_cd56bright)
de_markers_cd56bright$logP <- -log10(de_markers_cd56bright$p_val_adj + 1e-300)
de_markers_cd56bright$significant <- with(de_markers_cd56bright, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

top_up <- de_markers_cd56bright %>%
  filter(significant & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- de_markers_cd56bright %>%
  filter(significant & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

de_markers_cd56bright$regulation <- "Not Significant"
de_markers_cd56bright$regulation[de_markers_cd56bright$avg_log2FC >= 1 ] <- "Upregulated"
de_markers_cd56bright$regulation[de_markers_cd56bright$avg_log2FC <= -1] <- "Downregulated"

ggplot(de_markers_cd56bright, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = regulation), alpha = 0.6, size=2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not Significant" = "gray")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 5,
                  fontface = "bold",
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in CD56 Bright Subtype in CMV+ vs CMV- Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold"),  # Axis titles
        axis.text = element_text(size = 12, face = "bold")    # Axis tick labels
        
        )

####bright GSE analysis

# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- de_markers_cd56bright[!is.na(de_markers_cd56bright$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)
ranks <- -ranks
# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway), 
                                 levels = gsub("^HALLMARK_|^hallmark_", "", gsea_plot_data$pathway))

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:5]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:5]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(top_pathways$pathway, levels = top_pathways$pathway)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway), 
                               levels = gsub("^HALLMARK_|^hallmark_", "", top_pathways$pathway))

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))


# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))


ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: CD56 Bright in CMV+ v. CMV- Patients"
  ) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),     # Bold pathway labels
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold")
  )



######################  subtype CD56 proliferating ###################### 

cd56_proliferating <- subset(combined_seurat_subset, subset = annotation == "Proliferating")
Idents(cd56_proliferating) <- "cmv_status"

de_markers_proliferating <- FindMarkers(cd56_proliferating, ident.1 = "Positive", ident.2 = "Negative", 
                                     logfc.threshold = 0, min.pct = 0.1, assay = "RNA")

# Add gene names as a column
de_markers_proliferating$gene <- rownames(de_markers_proliferating)

# Classify significance and direction
de_markers_proliferating$significance <- "Not significant"
de_markers_proliferating$significance[de_markers_proliferating$avg_log2FC > 0.25 & de_markers_proliferating$p_val_adj < 0.05] <- "Upregulated" # 0.25 and 0.05
de_markers_proliferating$significance[de_markers_proliferating$avg_log2FC < -0.25 & de_markers_proliferating$p_val_adj < 0.05] <- "Downregulated"

# Factor levels to control legend order
de_markers_proliferating$significance <- factor(de_markers_proliferating$significance, levels = c("Upregulated", "Downregulated", "Not significant"))

de_markers_proliferating$gene <- rownames(de_markers_proliferating)
de_markers_proliferating$logP <- -log10(de_markers_proliferating$p_val_adj + 1e-300)
de_markers_proliferating$significant <- with(de_markers_proliferating, p_val_adj < 0.05 & abs(avg_log2FC) > 1) # 0.05 and 1

top_up <- de_markers_proliferating %>%
  filter(significant & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- de_markers_proliferating %>%
  filter(significant & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

de_markers_proliferating$regulation <- "Not Significant"
de_markers_proliferating$regulation[de_markers_proliferating$avg_log2FC >= 1 ] <- "Upregulated"
de_markers_proliferating$regulation[de_markers_proliferating$avg_log2FC <= -1] <- "Downregulated"

ggplot(de_markers_proliferating, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = regulation), alpha = 0.6, size=2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "Not Significant" = "gray")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in Proliferating Subtype in CMV+ vs CMV- Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")+
  theme(legend.position = "none")

####### proliferating GSEA ######


####bright GSE analysis

# now trying gsea in rstudio using the fgsea package
# Remove NA values and duplicates if any
markers_clean <- de_markers_proliferating[!is.na(de_markers_proliferating$avg_log2FC), ]
markers_clean <- markers_clean[!duplicated(rownames(markers_clean)), ]

# Create named vector: log2FC with gene names
ranks <- setNames(markers_clean$avg_log2FC, rownames(markers_clean))

# Sort in decreasing order
ranks <- sort(ranks, decreasing = TRUE)
ranks <- -ranks
# Install if needed
#BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

#BiocManager::install("fgsea", force=TRUE)
library(fgsea)

#BiocManager::install("BiocParallel", force = TRUE)

gsea_results <- fgsea(pathways = pathways,
                      stats = ranks,
                      minSize = 15,
                      maxSize = 500,
                      nperm = 10000
)

head(gsea_results[order(padj), ])

# Order results by NES
gsea_plot_data <- as.data.table(gsea_results)
gsea_plot_data <- gsea_plot_data[order(NES)]
gsea_plot_data$gs_name <- factor(gsea_plot_data$pathway, levels = gsea_plot_data$pathway)

# Select top 10 positively and negatively enriched pathways
top_up <- gsea_plot_data[order(-NES)][1:5]  # top 5 adaptive-enriched
top_down <- gsea_plot_data[order(NES)][1:5]  # top 5 CD56bright-enriched

# Combine the two
top_pathways <- rbind(top_up, top_down)

# Set factor levels to order the bars properly
top_pathways$gs_name <- factor(top_pathways$pathway, levels = top_pathways$pathway)

# Optional: Recalculate 'significant' column if needed
top_pathways$significant <- factor(top_pathways$padj < 0.05, levels = c(TRUE, FALSE))


# create significant definition
gsea_plot_data$significant <- factor(gsea_plot_data$padj < 0.05, levels = c(TRUE, FALSE))

ggplot(top_pathways, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(
    x = "Normalized Enrichment Score",
    y = "Pathway",
    title = "Top 5 Up & Downregulated Pathways: Proliferating NK in CMV+ v. CMV- Patients"
  )







################################

#######SUBSET VOLCANO PLOT#############


### SUBSET volcano plot !!!!!!

# set the Identity Class for Comparison
CMVpos_int$cmv_status <- "Positive"
CMVneg1$cmv_status <- "Negative"

# Merge into one Seurat object
combined_seurat <- merge(CMVpos_int, y = CMVneg1, add.cell.ids = c("CMVpos", "CMVneg1"))
Idents(combined_seurat) <- "cmv_status"

# Normalize the merged object

DefaultAssay(combined_seurat) <- "RNA"

combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- FindVariableFeatures(combined_seurat)

combined_seurat <- JoinLayers(combined_seurat)

markers_cmv <- FindMarkers(
  combined_seurat,
  ident.1 = "Positive",
  ident.2 = "Negative",
  logfc.threshold = 0,
  min.pct = 0.1,
  assay = "RNA"
)



markers_cmv$gene <- rownames(markers_cmv)
markers_cmv$logP <- -log10(markers_cmv$p_val_adj + 1e-300)
markers_cmv$significant <- with(markers_cmv, p_val_adj < 0.05 & abs(avg_log2FC) > 1)

library(ggplot2)
library(dplyr)

top_up <- markers_cmv %>%
  filter(significant & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_down <- markers_cmv %>%
  filter(significant & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

library(ggplot2)
library(ggrepel)  # For better text labels

library(ggplot2)
library(ggrepel)

ggplot(markers_cmv, aes(x = avg_log2FC, y = logP)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.5,
                  segment.color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: Gene Expression in CMV+ vs CMV− Patients",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 300)) +  # Adjust as needed
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")


# ok now let's plot those 5 upregulated genes on the umap and see where they fall
FeaturePlot(CMVpos, features = c("CD3D", "CD3E","MYOM2", "RPS26"))
Features(CMVpos)

# and downregulated genes - update genes here
FeaturePlot(CMVpos, features = c("IFI44L","MX1","XCL1", "IFITM3"))
Features(CMVpos)








############### MEMORY GSEA###################
# GSEA for memory-like phenotype

# Create a named vector of log2FC values, named by gene
ranked_genes <- markers_cmv %>%
  dplyr::select(gene, avg_log2FC) %>%
  arrange(desc(avg_log2FC)) %>%  # or arrange(avg_log2FC) if needed
  tibble::deframe()

# Check the output
head(ranked_genes)

# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Check the output
head(ranked_genes)

library(fgsea)
BiocManager::install("GSEABase", force=TRUE)
library(GSEABase)

# Load the gene set
gmt_file <- "ImmuneSigDB.gmt"
geneset <- getGmt(gmt_file)

# Convert gene set
pathways <- geneIds(geneset)

BiocManager::install("msigdbr", force=TRUE)
library(msigdbr)

# Example: Hallmark gene sets (for human)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# Run fgsea
fgsea_res <- fgsea(pathways = pathways, stats = ranked_genes, nperm = 1000)

# View results
head(fgsea_res)

# Order results by NES
gsea_plot_goldrath <- as.data.table(fgsea_res)
gsea_plot_goldrath <- gsea_plot_goldrath[order(NES)]
gsea_plot_goldrath$gs_name <- factor(gsea_plot_goldrath$pathway, levels = gsea_plot_goldrath$pathway)

# create significant definition
gsea_plot_goldrath$significant <- factor(gsea_plot_goldrath$padj < 0.05, levels = c(TRUE, FALSE))

gsea_plot_goldrath_significant <- gsea_plot_goldrath[gsea_plot_goldrath$significant == TRUE, ]

filtered_res <- gsea_plot_goldrath_significant %>%
  filter(str_detect(tolower(pathway), "memory|nk|pathway"))

# Optional: sort by NES or p-value if you want better plot ordering
filtered_res <- filtered_res %>% arrange(desc(NES))


ggplot(gsea_plot_goldrath, aes(x = NES, y = gs_name, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"), name = "padj < 0.05") +
  theme_minimal() +
  labs(x = "Normalized Enrichment Score", y = "Pathway",
       fill = "Significant", title = "Adaptive vs. CD56 Bright Pathway Enrichment")




################################

#export data

install.packages('devtools')
devtools::install_github('immunogenomics/presto')

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

# # Load necessary library
# library(openxlsx)
# 
# # Load necessary libraries
# library(Seurat)
# 
# # Assuming CMVpos is a list of two Seurat objects (e.g., CMVpos1 and CMVpos2)
# # Extract each Seurat object from the list:
# CMVpos1 <- CMVpos[1]
# CMVpos2 <- CMVpos[2]
# 
# # Function to extract data from each Seurat object
# extract_data_from_seurat <- function(seurat_obj) {
#   
#   variable_genes <- VariableFeatures(seurat_obj)
#   
#   # Extract RNA data (normalized expression data from the "integrated" assay)
#   rna_data <- GetAssayData(seurat_obj, assay = "integrated", slot = "data")
# 
#   # Extract PCA data (cell embeddings from PCA)
#   #pca_data <- seurat_obj[["pca"]]@cell.embeddings
#   
#   # Transpose to have cells as rows and genes as columns
#   rna_data <- as.data.frame(t(rna_data))
#   
#   # Extract UMAP data (cell embeddings from UMAP)
#   umap_data <- seurat_obj[["umap"]]@cell.embeddings
#   
#   # Extract cluster assignments
#   clusters <- seurat_obj$seurat_clusters
#   
#   # Extract metadata
#   metadata <- seurat_obj@meta.data
#   variable_gene_expression <- metadata[variable_genes, ]
#   
#   # Return all data as a list
#   return(list(
#     "RNA_Data" = as.data.frame(rna_data),
#     #"PCA_Data" = as.data.frame(pca_data),
#     "UMAP_Data" = as.data.frame(umap_data),
#     "Clusters" = data.frame(Cell_ID = names(clusters), Cluster = clusters),
#     "Metadata" = metadata,
#     "Variable_Gene_Expression" = as.data.frame(variable_gene_expression)
#   ))
# }
# 
# # Extract data from both Seurat objects (CMVpos1 and CMVpos2)
# data_CMVpos1 <- extract_data_from_seurat(CMVpos1)
# data_CMVpos2 <- extract_data_from_seurat(CMVpos2)
# 
# # Combine the data from both objects into a single list for export
# export_data <- list(
#   "CMVpos1_RNA_Data" = data_CMVpos1$RNA_Data,
#   "CMVpos2_RNA_Data" = data_CMVpos2$RNA_Data,
#   #"CMVpos1_PCA_Data" = data_CMVpos1$PCA_Data,
#   #"CMVpos2_PCA_Data" = data_CMVpos2$PCA_Data,
#   "CMVpos1_UMAP_Data" = data_CMVpos1$UMAP_Data,
#   "CMVpos2_UMAP_Data" = data_CMVpos2$UMAP_Data,
#   "CMVpos1_Clusters" = data_CMVpos1$Clusters,
#   "CMVpos2_Clusters" = data_CMVpos2$Clusters,
#   "CMVpos1_Metadata" = data_CMVpos1$Metadata,
#   "CMVpos2_Metadata" = data_CMVpos2$Metadata,
#   "CMVpos1_Var_Gene_Expression" = data_CMVpos1$Variable_Gene_Expression,
#   "CMVpos2_Var_Gene_Expression" = data_CMVpos2$Variable_Gene_Expression
# )
# 
# # Write the combined data to an Excel file
# write.xlsx(export_data, file = "CMVpos_integrated_data.xlsx")

# 
# # Name populations based on hashtags
# CMVpos$population <- "NKG2Cneg"
# CMVpos$population[CMVpos$HTO_maxID=="NKG2Cpos"] <- "NKG2Cpos"
# CMVpos$population[CMVpos$HTO_corrected_maxID%in%c("CMVpos1-NKG2Cpos","CMVpos5-NKG2Cpos")] <- "NKG2Cpos"
# 
# DimPlot(CMVpos, group.by = "population", cols = c("grey",colorscale[7]), shuffle = TRUE)&NoLegend()
# 
# # Integration of CMVneg (use the same anchor features as for CMVpos to enable merging afterwards)
# integration_anchors_CMVneg <- FindIntegrationAnchors(c(my_SeuratList[[4]], my_SeuratList[[5]]),
#                                                      dims = 1:10, anchor.features = shared_variable_features_CMVpos)
# CMVneg <- IntegrateData(integration_anchors_CMVneg, dims = 1:10, features.to.integrate = all_features)

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

library(viridisLite)  # or library(viridis)
viridis_colors <- viridis(2)  # or more, e.g., viridis(10) for smoother gradient

# gene expression
FeaturePlot(CMVneg, features = c("CX3CR1","GZMK","STMN1"), cols = viridis_colors,
            min.cutoff=0,
            max.cutoff=4)
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






