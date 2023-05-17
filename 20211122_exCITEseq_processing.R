library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(plotly)
library(Seurat)

setwd("~/OneDrive - NYU Langone Health/ECCITE-seq/Experiments/20211114_ECCITEseq Pilot_Subject2")

#------------------------------------Initial QC-----------------------------------

stim.out <- Read10X(data.dir = "Stim_outs_raw")
unstim.out <- Read10X(data.dir = "Unstim_outs_raw")

#creating Seurat object
stim <- CreateSeuratObject(counts = stim.out$`Gene Expression`,  Project = "stim")
unstim <- CreateSeuratObject(counts = unstim.out$`Gene Expression`,  Project = "unstim")

#new assay to store ADT information
adt_stim <- stim.out$`Antibody Capture`
adt_unstim <- unstim.out$`Antibody Capture`

#remove ADTs that were not put into the panel 
adt_assay_stim <- adt_stim[setdiff(rownames(x = adt_stim),c("7moPostVax", "1moPostBoost", "CD69.1", "CD137_4-1BB", "CD154")), ]
adt_assay_unstim <- adt_unstim[setdiff(rownames(x = adt_unstim),c("7moPostVax", "1moPostBoost", "CD69.1","CD137_4-1BB", "CD154")), ]

adt_assay_stim <- CreateAssayObject(counts = adt_assay_stim)
adt_assay_unstim <- CreateAssayObject(counts = adt_assay_unstim)

#add assay to previously created Seurat object
stim[["ADT"]] <- adt_assay_stim
tetramer[["ADT"]] <- adt_assay_tetramer
whole[["ADT"]] <- adt_assay_whole
unstim[["ADT"]] <- adt_assay_unstim

#add HTO as separate assay
HTO_assay_stim <- adt_stim[intersect(rownames(x = adt_stim),c("7moPostVax", "1moPostBoost")), ]
HTO_assay_unstim <- adt_unstim[intersect(rownames(x = adt_unstim),c("7moPostVax", "1moPostBoost")), ]

HTO_assay_stim <- CreateAssayObject(counts = HTO_assay_stim)
HTO_assay_unstim <- CreateAssayObject(counts = HTO_assay_unstim)

#add assay to previously created Seurat object
stim[["HTO"]] <- HTO_assay_stim
unstim[["HTO"]] <- HTO_assay_unstim

#validate object has multiple assays
Assays(stim)
rownames(stim[["HTO"]])
rownames(stim[["ADT"]])

#set identity for each lane 
stim@meta.data$orig.ident <- "stim"
unstim@meta.data$orig.ident <- "unstim"

#Barcode rank plots across three modalities
stim_Barcode_rna <- CalculateBarcodeInflections(stim, group.column = 'orig.ident', barcode.column = "nCount_RNA")
s_rna <- BarcodeInflectionsPlot(stim_Barcode_rna) + scale_x_log10() +scale_y_log10()
stim_Barcode_adt <- CalculateBarcodeInflections(stim, group.column = 'orig.ident', barcode.column = "nCount_ADT")
s_adt <- BarcodeInflectionsPlot(stim_Barcode_adt) + scale_x_log10() +scale_y_log10() 
stim_Barcode_hto <- CalculateBarcodeInflections(stim, group.column = 'orig.ident', barcode.column = "nCount_HTO")
s_hto <- BarcodeInflectionsPlot(stim_Barcode_hto) + scale_x_log10() +scale_y_log10()

s_rna
s_adt
s_hto


unstim_Barcode_rna <- CalculateBarcodeInflections(unstim, group.column = 'orig.ident', barcode.column = "nCount_RNA")
u_rna <- BarcodeInflectionsPlot(unstim_Barcode_rna) + scale_x_log10() +scale_y_log10() 
unstim_Barcode_adt <- CalculateBarcodeInflections(unstim, group.column = 'orig.ident', barcode.column = "nCount_ADT")
u_adt <- BarcodeInflectionsPlot(unstim_Barcode_adt) + scale_x_log10() +scale_y_log10() 
unstim_Barcode_hto <- CalculateBarcodeInflections(unstim, group.column = 'orig.ident', barcode.column = "nCount_HTO")
u_hto <- BarcodeInflectionsPlot(unstim_Barcode_hto) + scale_x_log10() +scale_y_log10()

u_rna
u_adt
u_hto

#separate out empty droplets based on low nCount 
empty_droplets_stim <- subset(x = stim, subset = nCount_RNA <= 1500)
empty_droplets_unstim <- subset(x = unstim, subset = nCount_RNA <= 1500)


empty_droplets_stim <- Matrix::rowSums(empty_droplets_stim@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim$marker<- rownames(empty_droplets_stim)
names(empty_droplets_stim)[1] <- 'Counts'
empty_droplets_stim$sample<- "empty_droplets_stim"
empty_droplets_stim$lane <- "stim_empty"
empty_droplets_stim$type <- "empty"

empty_droplets_unstim <- Matrix::rowSums(empty_droplets_unstim@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_unstim$marker<- rownames(empty_droplets_unstim)
names(empty_droplets_unstim)[1] <- 'Counts'
empty_droplets_unstim$sample<- "empty_droplets_unstim"
empty_droplets_unstim$lane <- "unstim_empty"
empty_droplets_unstim$type <- "empty"


#separate out cells  based on high nCount 
cells_stim <- subset(x = stim, subset = nCount_RNA > 1500)
cells_unstim <- subset(x = unstim, subset = nCount_RNA > 1500)

cells_stim <- Matrix::rowSums(cells_stim@assays[["ADT"]]) %>% as.data.frame() 
cells_stim$marker<- rownames(cells_stim)
names(cells_stim)[1] <- 'Counts'
cells_stim$sample<- "cells_stim"
cells_stim$lane <- "stim"
cells_stim$type <- cells_stim$marker

cells_unstim <- Matrix::rowSums(cells_unstim@assays[["ADT"]]) %>% as.data.frame() 
cells_unstim$marker<- rownames(cells_unstim)
names(cells_unstim)[1] <- 'Counts'
cells_unstim$sample<- "cells_unstim"
cells_unstim$lane <- "unstim"
cells_unstim$type <- cells_unstim$marker

all_stim <- rbind(empty_droplets_stim,  cells_stim)
all_unstim <- rbind(empty_droplets_unstim, cells_unstim)

Marker_order <- empty_droplets_stim  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim <- Marker_order$marker

Marker_order <- empty_droplets_unstim  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_unstim <- Marker_order$marker

all_stim$marker <- factor(all_stim$marker, levels = order_stim)
all_unstim$marker <- factor(all_unstim$marker, levels = order_unstim)

colourCount = length(unique(all_stim$type))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

Sti <- ggplot(all_stim, aes( x = "UMI", y = Counts/10^6, fill=marker)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") + facet_wrap(~lane) + 
  scale_fill_manual(values = getPalette(colourCount))

Uns <- ggplot(all_unstim, aes( x = "UMI", y = Counts/10^6, fill=marker)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") + facet_wrap(~lane) + 
  scale_fill_manual(values = getPalette(colourCount))

Sti  + Uns


Sti <- ggplot(all_stim, aes( x = marker, y = Counts/10^6, fill = sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("darkorchid",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim") + theme(legend.position="top")

Uns <- ggplot(all_unstim, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("blue",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Unstim") + theme(legend.position="top")

Sti + coord_flip() + labs(title= "Stim")
Uns + coord_flip() + labs(title= "Unstim")


Sti + coord_flip() | Uns + coord_flip()

#-----------------------------------------objects moving forward-----------------------------------------

cells_stim <- subset(x = stim, subset = nCount_RNA > 500 & nFeature_RNA > 400)
cells_unstim <- subset(x = unstim, subset = nCount_RNA > 500 & nFeature_RNA > 400)

#-----------------------------------------demultiplex?--------------------------------------------------
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cells_stim <- NormalizeData(cells_stim, assay = "HTO", normalization.method = "CLR")
cells_unstim <- NormalizeData(cells_unstim, assay = "HTO", normalization.method = "CLR")

#We perform a k-medoid clustering on the normalized HTO values, which initially separates cells into K(# of samples)+1 clusters.
#We calculate a ‘negative’ distribution for HTO. For each HTO, we use the cluster with the lowest average value as the negative group.
#For each HTO, we fit a negative binomial distribution to the negative cluster. We use the 0.99 quantile of this distribution as a threshold.
#Based on these thresholds, each cell is classified as positive or negative for each HTO.
#Cells that are positive for more than one HTOs are annotated as doublets.

#positive.quantile - The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#Determine how less stringent demultiplex would look by the highest HTO signal at or above a set UMI threshold
cells_stim <- HTODemux(cells_stim, assay = "HTO", positive.quantile = 0.99)
table(cells_stim$HTO_classification.global) 
table(cells_stim$HTO_maxID) 
# Group cells based on the max HTO signal
Idents(cells_stim) <- "HTO_maxID"
RidgePlot(cells_stim, assay = "HTO", features = rownames(cells_stim[["HTO"]])[1:2], ncol = 2)
s <- FeatureScatter(cells_stim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="orig.ident") + labs(title = "stim") + NoLegend()
Idents(cells_stim) <- "HTO_classification.global"
VlnPlot(cells_stim, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
RidgePlot(cells_stim, assay = "HTO", features = rownames(cells_stim[["HTO"]]), group.by="HTO_classification.global", combine=TRUE)

cells_unstim <- HTODemux(cells_unstim, assay = "HTO", positive.quantile = 0.99)
table(cells_unstim$HTO_classification.global)
# Group cells based on the max HTO signal
Idents(cells_unstim) <- "HTO_maxID"
RidgePlot(cells_unstim, assay = "HTO", features = rownames(cells_unstim[["HTO"]])[1:2], ncol = 2)
u <- FeatureScatter(cells_unstim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost",group.by="orig.ident") + labs(title = "unstim")
Idents(cells_unstim) <- "HTO_classification.global"
VlnPlot(cells_unstim, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#FEATURE SCATTER PLOT OF HTOS 
s| u 

## Determine how less stringent demultiplex would look with HTO_maxID at or above a set UMI threshold of 0.25
cells_stim <- subset(cells_stim, subset = hto_7moPostVax > 0.25 | hto_1moPostBoost > 0.25)
cells_unstim <- subset(cells_unstim, subset = hto_7moPostVax > 0.25 | hto_1moPostBoost > 0.25)

stim_ID <- cells_stim$HTO_classification.global
stim_ID[cells_stim@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim$stim_ID <- factor(stim_ID)

unstim_ID <- cells_unstim$HTO_classification.global
unstim_ID[cells_unstim@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_unstim$unstim_ID <- factor(unstim_ID)

s <- FeatureScatter(cells_stim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="stim_ID") + labs(title = "stim") + NoLegend()
u <- FeatureScatter(cells_unstim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost",group.by="unstim_ID") + labs(title = "unstim")

RidgePlot(cells_tetramer, assay = "HTO", features = rownames(cells_stim[["HTO"]]), group.by="hash.ID", combine=TRUE)

#FEATURE SCATTER PLOT OF HTOS with less stringent demultiplex 
s | u 

#look at numbers of singlets and doublets
stim_hash <- cells_stim@meta.data %>% 
  ggplot(aes(x=stim_ID)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim NCells") #+ ylim(c(0,30000))


unstim_hash <-cells_unstim@meta.data %>% 
  ggplot(aes(x=unstim_ID)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells")# + ylim(c(0,5000))

stim_hash  | unstim_hash

#-----------------------------------------filter out doublets--------------------------------------------------
library(scDblFinder)
DefaultAssay(cells_stim) <- "RNA"
DefaultAssay(cells_unstim) <- "RNA"

sce.stim <- as.SingleCellExperiment(cells_stim)
sce.stim <- scDblFinder(sce.stim, samples="HTO_maxID")
table(sce.stim$scDblFinder.class, sce.stim$HTO_maxID)

sce.uns <- as.SingleCellExperiment(cells_unstim)
sce.uns <- scDblFinder(sce.uns, samples="HTO_maxID")
table(sce.uns$scDblFinder.class, sce.uns$HTO_maxID)

identical(colnames(cells_stim),colnames(sce.stim))
identical(colnames(cells_unstim),colnames(sce.uns))

cells_stim$scDblFinder.class <- sce.stim$scDblFinder.class
cells_stim$scDblFinder.score <- sce.stim$scDblFinder.score
cells_stim$scDblFinder.weighted <- sce.stim$scDblFinder.weighted
rm(sce.stim)

cells_unstim$scDblFinder.class <- sce.uns$scDblFinder.class
cells_unstim$scDblFinder.score <- sce.uns$scDblFinder.score
cells_unstim$scDblFinder.weighted <- sce.uns$scDblFinder.weighted
rm(sce.uns)

cells_stim$scDblFinder.class <- factor(cells_stim$scDblFinder.class, levels=c("singlet","doublet"))
cells_unstim$scDblFinder.class <- factor(cells_unstim$scDblFinder.class, levels=c("singlet","doublet"))

s <- FeatureScatter(cells_stim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by ="scDblFinder.class") + labs(title = "stim")
u <- FeatureScatter(cells_unstim, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="scDblFinder.class") + labs(title = "unstim")

s | u

cells_stim$HTO_maxID <- factor(cells_stim$HTO_maxID, levels=c("7moPostVax","1moPostBoost"))
cells_unstim$HTO_maxID <- factor(cells_unstim$HTO_maxID, levels=c("7moPostVax","1moPostBoost"))


#RENAME HTO TO Subject2_7moPostVax and Subject2_1moPostBoost

#-----------------------------------------REMOVE UNUSED OBJECTS-----------------------------------------
remove(stim.out,tetramer.out,whole.out, unstim.out, adt_stim,adt_tetramer,adt_whole,adt_unstim,stim, tetramer, whole, unstim, 
       stim_Barcode_adt, stim_Barcode_hto, stim_Barcode_rna, tet_Barcode_adt, tet_Barcode_hto, tet_Barcode_rna,
       who_Barcode_adt, who_Barcode_hto, who_Barcode_rna, unstim_Barcode_adt, unstim_Barcode_hto, unstim_Barcode_rna, 
       adt_assay_stim, adt_assay_tetramer, adt_assay_whole, HTO_assay_stim, HTO_assay_tetramer, HTO_assay_whole, HTO_assay_unstim)

#-----------------------------------------stim-----------------------------------------
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim <- AddMetaData(object = cells_stim, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) + labs(title = "stim") 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim, feature1 = "nCount_RNA", feature2 = "percent.mito") + geom_hline(yintercept = 0.08) + labs(title = "stim")
FeatureScatter(object = cells_stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "scDblFinder.class") +  geom_hline(yintercept = c(400, 6000)) + labs(title = "stim")

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim <- subset(x = cells_stim, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim <- NormalizeData(cells_stim, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim <- FindVariableFeatures(cells_stim, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim <- ScaleData(cells_stim, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim <- RunPCA(cells_stim, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim, reduction = "pca")
DimHeatmap(object = cells_stim, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)

cells_stim <- FindNeighbors(cells_stim, dims = 1:30)
cells_stim <- FindClusters(cells_stim, resolution = 0.8, verbose = FALSE)
cells_stim <- RunUMAP(cells_stim, dims = 1:30)
DimPlot(cells_stim, label = TRUE, group.by = "scDblFinder.class", pt.size = 1)

# Normalize ADT data,
cells_stim <- NormalizeData(cells_stim, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim <- ScaleData(cells_stim, assay = "ADT")


#-----------------------------------------unstim-----------------------------------------
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_unstim@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_unstim@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_unstim@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_unstim <- AddMetaData(object = cells_unstim, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_unstim, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident") + labs(title = "unstim")

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_unstim, feature1 = "nCount_RNA", feature2 = "percent.mito") + geom_hline(yintercept = 0.08) + labs(title = "unstim")
FeatureScatter(object = cells_unstim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "scDblFinder.class") +  geom_hline(yintercept = c(200, 3500)) + labs(title = "unstim")

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_unstim <- subset(x = cells_unstim, subset = nFeature_RNA > 400 & nFeature_RNA < 3500 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_unstim <- NormalizeData(cells_unstim, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_unstim <- FindVariableFeatures(cells_unstim, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_unstim))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_unstim <- ScaleData(cells_unstim, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_unstim <- RunPCA(cells_unstim, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_unstim, reduction = "pca")
DimHeatmap(object = cells_unstim, reduction = "pca", cells = 500, balanced = TRUE)

cells_unstim <- FindNeighbors(cells_unstim, dims = 1:30)
cells_unstim <- FindClusters(cells_unstim, resolution = 0.8, verbose = FALSE)
cells_unstim <- RunUMAP(cells_unstim, dims = 1:30)
DimPlot(cells_unstim, label = FALSE, group.by = "scDblFinder.class")

# Normalize ADT data,
cells_unstim <- NormalizeData(cells_unstim, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_unstim <- ScaleData(cells_unstim, assay = "ADT")




#-----------------------------------Filtering out doublets?-------------------------------------------
#look at numbers of signlets and doublets
stim_hash <- cells_stim@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim NCells") + ylim(c(0,30000))

unstim_hash <-cells_unstim@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells") + ylim(c(0,5000))
  
stim_hash  | unstim_hash

#CD19 vs CD3 gene expression before filtering out doublets
stim_hash <- FeatureScatter(cells_stim, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") #+ NoLegend() #+ geom_vline(xintercept = 0.75) +  geom_hline(yintercept = 0.75)
unstim_hash <- FeatureScatter(cells_unstim, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class")# + NoLegend() #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)

stim_hash| unstim_hash

#subset to on singlets 
stim_singlet <- subset(cells_stim, subset=scDblFinder.class=="singlet")
uns_singlet <- subset(cells_unstim, subset=scDblFinder.class=="singlet")

s <- FeatureScatter(stim_singlet, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by ="HTO_maxID") + labs(title = "stim")
u <- FeatureScatter(uns_singlet, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="HTO_maxID") + labs(title = "unstim")

s | u

#see what doublets look like 
stim_doub <- subset(cells_stim, subset=scDblFinder.class=="doublet")
uns_doub <- subset(cells_unstim, subset=scDblFinder.class=="doublet")
 
s <- FeatureScatter(stim_doub, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by ="HTO_maxID") + labs(title = "stim_doublet")
u <- FeatureScatter(uns_doub, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="HTO_maxID") + labs(title = "unstim_doublet")

s | u

#CD19 vs CD3 gene expression for doublets
stim_hash <- FeatureScatter(stim_doub, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class")+ labs(title = "stim_doublet") #+ NoLegend() #+ geom_vline(xintercept = 0.75) +  geom_hline(yintercept = 0.75)
unstim_hash <- FeatureScatter(uns_doub, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + labs(title = "unstim_doublet")# + NoLegend() #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)

stim_hash | unstim_hash

#subset singlets round 2 by excluding cells that are double positive for CD19 and CD3E genes 
stim_hash <- FeatureScatter(stim_singlet, feature1 = "CD3E", feature2 = "CD19", group.by = "HTO_maxID") + geom_vline(xintercept = 0.25) +  geom_hline(yintercept = 0.25) + labs(title = "stim")
unstim_hash <- FeatureScatter(uns_singlet, feature1 = "CD3E", feature2 = "CD19", group.by = "HTO_maxID")+geom_vline(xintercept = 0.25)+ geom_hline(yintercept = 0.25) + labs(title = "unstim")

stim_hash | unstim_hash

stim_singlet.2 <- subset(stim_singlet, subset = CD3E < 0.25 & CD19 > 0.25 |CD3E > 0.25 & CD19 < 0.25 | CD3E < 0.25 & CD19 < 0.25)
uns_singlet.2 <- subset(uns_singlet, subset = CD3E < 0.25 & CD19 > 0.25 |CD3E > 0.25 & CD19 < 0.25| CD3E < 0.25 & CD19 < 0.25 )

stim_hash <- FeatureScatter(stim_singlet.2, feature1 = "CD3E", feature2 = "CD19", group.by = "HTO_maxID") + geom_vline(xintercept = 0.25) +  geom_hline(yintercept = 0.25) + labs(title = "stim")
unstim_hash <- FeatureScatter(uns_singlet.2, feature1 = "CD3E", feature2 = "CD19", group.by = "HTO_maxID")+geom_vline(xintercept = 0.25)+ geom_hline(yintercept = 0.25) + labs(title = "unstim")

stim_hash  | unstim_hash

#see if cleans up demultiplexing
s <- FeatureScatter(stim_singlet.2, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by ="HTO_maxID") + labs(title = "stim")
u <- FeatureScatter(uns_singlet.2, feature1 = "hto_7moPostVax", feature2 = "hto_1moPostBoost", group.by="HTO_maxID") + labs(title = "unstim")

s | u

#look at number of singlets
stim_hash <- stim_singlet.2@meta.data %>% 
  ggplot(aes(x=HTO_maxID)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim NCells")

unstim_hash <-uns_singlet.2@meta.data %>% 
  ggplot(aes(x=HTO_maxID)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells")

stim_hash | unstim_hash

#------------------------------------SAVING final objects--------------------------------------------------
stim_singlet.2@meta.data$subject <- "Subject2"
stim_singlet.2@meta.data$timepoint <- stim_singlet.2@meta.data$HTO_maxID
stim_singlet.2@meta.data <- stim_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject2") ~ "no"
  ))
stim_singlet.2$timepoint <- factor(stim_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost"))
View(stim_singlet.2@meta.data)

uns_singlet.2@meta.data$subject <- "Subject2"
uns_singlet.2@meta.data$timepoint <- uns_singlet.2@meta.data$HTO_maxID
uns_singlet.2@meta.data <- uns_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject2") ~ "no"
  ))
uns_singlet.2$timepoint <- factor(uns_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost"))

setwd("~/OneDrive - NYU Langone Health/ECCITE-seq/SavedObjects")
saveRDS(stim_singlet.2, file = "20211114_stim.rds")
saveRDS(uns_singlet.2, file = "20211114_unstim.rds")
