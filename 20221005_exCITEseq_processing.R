library(ggplot2)
library(tidyverse)
library(stringr)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(gplots)
library(plotly)
library(Seurat)
sessionInfo()

setwd("/Volumes/heratr01lab/heratr01labspace/ECCITE-seq/Differential_Spike_Priming_Studies/Experiments/20220927_ECCITEseq_AIM")
------------------------------------Initial QC-----------------------------------

stim1.out <- Read10X(data.dir = "Stim1_outs_raw")
stim2.out <- Read10X(data.dir = "Stim2_outs_raw")
stim3.out <- Read10X(data.dir = "Stim3_outs_raw")
stim4.out <- Read10X(data.dir = "Stim4_outs_raw")
stim5.out <- Read10X(data.dir = "Stim5_outs_raw")
unstim.out <- Read10X(data.dir = "Unstim_outs_raw")

#creating Seurat object
stim1 <- CreateSeuratObject(counts = stim1.out$`Gene Expression`,  Project = "stim1")
stim2 <- CreateSeuratObject(counts = stim2.out$`Gene Expression`,  Project = "stim2")
stim3 <- CreateSeuratObject(counts = stim3.out$`Gene Expression`,  Project = "stim3")
stim4 <- CreateSeuratObject(counts = stim4.out$`Gene Expression`,  Project = "stim4")
stim5 <- CreateSeuratObject(counts = stim5.out$`Gene Expression`,  Project = "stim5")
unstim <- CreateSeuratObject(counts = unstim.out$`Gene Expression`,  Project = "unstim")

#new assay to store ADT information
adt_stim1 <- stim1.out$`Antibody Capture`
adt_stim2 <-  stim2.out$`Antibody Capture`
adt_stim3 <- stim3.out$`Antibody Capture`
adt_stim4 <- stim4.out$`Antibody Capture`
adt_stim5 <- stim5.out$`Antibody Capture`
adt_unstim <- unstim.out$`Antibody Capture`
adt_stim1@Dimnames[[1]]
#remove HTOs from ADT panel 
adt_assay_stim1 <- adt_stim1[setdiff(rownames(x = adt_stim1),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
adt_assay_stim2 <- adt_stim2[setdiff(rownames(x = adt_stim2),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
adt_assay_stim3 <- adt_stim3[setdiff(rownames(x = adt_stim3),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
adt_assay_stim4 <- adt_stim4[setdiff(rownames(x = adt_stim4),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
adt_assay_stim5 <- adt_stim5[setdiff(rownames(x = adt_stim5),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
adt_assay_unstim <- adt_unstim[setdiff(rownames(x = adt_unstim),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                   "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                   "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]

adt_assay_stim1 <- CreateAssayObject(counts = adt_assay_stim1)
adt_assay_stim2 <- CreateAssayObject(counts = adt_assay_stim2)
adt_assay_stim3 <- CreateAssayObject(counts = adt_assay_stim3)
adt_assay_stim4 <- CreateAssayObject(counts = adt_assay_stim4)
adt_assay_stim5 <- CreateAssayObject(counts = adt_assay_stim5)
adt_assay_unstim <- CreateAssayObject(counts = adt_assay_unstim)

#add assay to previously created Seurat object
stim1[["ADT"]] <- adt_assay_stim1
stim2[["ADT"]] <- adt_assay_stim2
stim3[["ADT"]] <- adt_assay_stim3
stim4[["ADT"]] <- adt_assay_stim4
stim5[["ADT"]] <- adt_assay_stim5
unstim[["ADT"]] <- adt_assay_unstim

#add HTO as separate assay
HTO_assay_stim1 <- adt_stim1[intersect(rownames(x = adt_stim1),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                  "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                  "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
HTO_assay_stim2 <- adt_stim2[intersect(rownames(x = adt_stim2),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                  "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                  "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
HTO_assay_stim3 <- adt_stim3[intersect(rownames(x = adt_stim3),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                  "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                  "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
HTO_assay_stim4 <- adt_stim4[intersect(rownames(x = adt_stim4),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                  "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                  "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
HTO_assay_stim5 <- adt_stim5[intersect(rownames(x = adt_stim5),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru","Subject4_PostBreakThru",
                                                                  "Subject4_7moPostVax","Subject4_1moPostBoost","Subject10_7moPostVax","Subject10_1moPostBoost",
                                                                  "Subject11_7moPostVax","Subject11_1moPostBoost","Subject5_PostBreakThru")), ]
HTO_assay_unstim <- adt_unstim[intersect(rownames(x = adt_unstim),c( "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru",
                                                                     "Subject10_1moPostBoost","Subject11_7moPostVax","Subject11_1moPostBoost")), ]

HTO_assay_stim1 <- CreateAssayObject(counts = HTO_assay_stim1)
HTO_assay_stim2 <- CreateAssayObject(counts = HTO_assay_stim2)
HTO_assay_stim3 <- CreateAssayObject(counts = HTO_assay_stim3)
HTO_assay_stim4 <- CreateAssayObject(counts = HTO_assay_stim4)
HTO_assay_stim5 <- CreateAssayObject(counts = HTO_assay_stim5)
HTO_assay_unstim <- CreateAssayObject(counts = HTO_assay_unstim)

#add assay to previously created Seurat object
stim1[["HTO"]] <- HTO_assay_stim1
stim2[["HTO"]] <- HTO_assay_stim2
stim3[["HTO"]] <- HTO_assay_stim3
stim4[["HTO"]] <- HTO_assay_stim4
stim5[["HTO"]] <- HTO_assay_stim5
unstim[["HTO"]] <- HTO_assay_unstim

#validate object has multiple assays
Assays(stim1)
rownames(stim1[["HTO"]])
rownames(unstim[["HTO"]])

#set identity for each lane 
stim1@meta.data$orig.ident <- "stim1"
stim2@meta.data$orig.ident <- "stim2"
stim3@meta.data$orig.ident <- "stim3"
stim4@meta.data$orig.ident <- "stim4"
stim5@meta.data$orig.ident <- "stim5"
unstim@meta.data$orig.ident <- "unstim"

saveRDS(stim1, file="stim1_raw.rds")
saveRDS(stim2, file="stim2_raw.rds")
saveRDS(stim3, file="stim3_raw.rds")
saveRDS(stim4, file="stim4_raw.rds")
saveRDS(stim5, file="stim5_raw.rds")
saveRDS(unstim, file="unstim_raw.rds")


#Barcode rank plots across three modalities
stim1 <- CalculateBarcodeInflections(stim1, group.column = 'orig.ident', barcode.column = "nCount_RNA", threshold.low = 48000, threshold.high = 53000)
s1_rna <- BarcodeInflectionsPlot(stim1) + scale_x_log10() +scale_y_log10()
s1_rna

stim2 <- CalculateBarcodeInflections(stim2, group.column = 'orig.ident', barcode.column = "nCount_RNA", threshold.low = 40000, threshold.high = 45000)
s2_rna <- BarcodeInflectionsPlot(stim2) + scale_x_log10() +scale_y_log10()
s2_rna

stim3 <- CalculateBarcodeInflections(stim3, group.column = 'orig.ident', barcode.column = "nCount_RNA", threshold.low = 24000, threshold.high = 29000)
s3_rna <- BarcodeInflectionsPlot(stim3) + scale_x_log10() +scale_y_log10()
s3_rna

stim4 <- CalculateBarcodeInflections(stim4, group.column = 'orig.ident', barcode.column = "nCount_RNA", threshold.low = 38000, threshold.high = 43000)
s4_rna <- BarcodeInflectionsPlot(stim4) + scale_x_log10() +scale_y_log10()
s4_rna

stim5 <- CalculateBarcodeInflections(stim5, group.column = 'orig.ident', barcode.column = "nCount_RNA", threshold.low = 32000, threshold.high = 37000)
s5_rna <- BarcodeInflectionsPlot(stim5) + scale_x_log10() +scale_y_log10()
s5_rna

unstim <- CalculateBarcodeInflections(unstim, group.column = 'orig.ident', barcode.column = "nCount_RNA",threshold.low = 23000 , threshold.high = 27000)
u_rna <- BarcodeInflectionsPlot(unstim)  + scale_x_log10() +scale_y_log10()
u_rna

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
Rank_s1 <- stim1@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
stim1 <- AddMetaData(object = stim1, metadata = Rank_s1, col.name = c( "nCount_RNA.1","rank"))

Rank_s2 <- stim2@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
stim2 <- AddMetaData(object = stim2, metadata = Rank_s2, col.name = c( "nCount_RNA.1","rank"))

Rank_s3 <- stim3@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
stim3 <- AddMetaData(object = stim3, metadata = Rank_s3, col.name = c( "nCount_RNA.1","rank"))

Rank_s4 <- stim4@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
stim4 <- AddMetaData(object = stim4, metadata = Rank_s4, col.name = c( "nCount_RNA.1","rank"))

Rank_s5 <- stim5@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
stim5 <- AddMetaData(object = stim5, metadata = Rank_s5, col.name = c( "nCount_RNA.1","rank"))

Rank_u <- unstim@tools$CalculateBarcodeInflections$barcode_distribution[,2:3]
unstim <- AddMetaData(object = unstim, metadata = Rank_u, col.name = c( "nCount_RNA.1","rank"))

#separate out empty droplets based on inflection point determined in CalculateBarcodeInflections
empty_droplets_stim1 <- subset(x = stim1, subset = rank > 2*stim1@tools$CalculateBarcodeInflections$inflection_points$rank)
empty_droplets_stim2 <- subset(x = stim2, subset = rank > 2*stim2@tools$CalculateBarcodeInflections$inflection_points$rank)
empty_droplets_stim3 <- subset(x = stim3, subset = rank > 2*stim3@tools$CalculateBarcodeInflections$inflection_points$rank)
empty_droplets_unstim <- subset(x = unstim, subset = rank > 2*unstim@tools$CalculateBarcodeInflections$inflection_points$rank)

empty_droplets_stim1 <- Matrix::rowSums(empty_droplets_stim1@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim1$marker<- rownames(empty_droplets_stim1)
names(empty_droplets_stim1)[1] <- 'Counts'
empty_droplets_stim1$sample<- "empty_droplets_stim1"
empty_droplets_stim1$lane <- "stim1_empty"
empty_droplets_stim1$type <- "empty"

empty_droplets_stim2 <- Matrix::rowSums(empty_droplets_stim2@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim2$marker<- rownames(empty_droplets_stim2)
names(empty_droplets_stim2)[1] <- 'Counts'
empty_droplets_stim2$sample<- "empty_droplets_stim2"
empty_droplets_stim2$lane <- "stim2_empty"
empty_droplets_stim2$type <- "empty"

empty_droplets_stim3 <- Matrix::rowSums(empty_droplets_stim3@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim3$marker<- rownames(empty_droplets_stim3)
names(empty_droplets_stim3)[1] <- 'Counts'
empty_droplets_stim3$sample<- "empty_droplets_stim3"
empty_droplets_stim3$lane <- "stim3_empty"
empty_droplets_stim3$type <- "empty"

empty_droplets_stim4 <- Matrix::rowSums(empty_droplets_stim4@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim4$marker<- rownames(empty_droplets_stim4)
names(empty_droplets_stim4)[1] <- 'Counts'
empty_droplets_stim4$sample<- "empty_droplets_stim4"
empty_droplets_stim4$lane <- "stim4_empty"
empty_droplets_stim4$type <- "empty"

empty_droplets_stim5 <- Matrix::rowSums(empty_droplets_stim5@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_stim5$marker<- rownames(empty_droplets_stim5)
names(empty_droplets_stim5)[1] <- 'Counts'
empty_droplets_stim5$sample<- "empty_droplets_stim5"
empty_droplets_stim5$lane <- "stim5_empty"
empty_droplets_stim5$type <- "empty"

empty_droplets_unstim <- Matrix::rowSums(empty_droplets_unstim@assays[["ADT"]]) %>% as.data.frame() 
empty_droplets_unstim$marker<- rownames(empty_droplets_unstim)
names(empty_droplets_unstim)[1] <- 'Counts'
empty_droplets_unstim$sample<- "empty_droplets_unstim"
empty_droplets_unstim$lane <- "unstim_empty"
empty_droplets_unstim$type <- "empty"

#separate out cells  based on high nCount 
cells_stim1 <- subset(x = stim1, subset = rank <= stim1@tools$CalculateBarcodeInflections$inflection_points$rank)
cells_stim2 <- subset(x = stim2, subset = rank <= stim2@tools$CalculateBarcodeInflections$inflection_points$rank)
cells_stim3 <- subset(x = stim3, subset = rank <= stim3@tools$CalculateBarcodeInflections$inflection_points$rank)
cells_stim4 <- subset(x = stim4, subset = rank <= stim4@tools$CalculateBarcodeInflections$inflection_points$rank)
cells_stim5 <- subset(x = stim5, subset = rank <= stim5@tools$CalculateBarcodeInflections$inflection_points$rank)
cells_unstim <- subset(x = unstim, subset = rank <= unstim@tools$CalculateBarcodeInflections$inflection_points$rank)

cells_stim1 <- Matrix::rowSums(cells_stim1@assays[["ADT"]]) %>% as.data.frame() 
cells_stim1$marker<- rownames(cells_stim1)
names(cells_stim1)[1] <- 'Counts'
cells_stim1$sample<- "cells_stim1"
cells_stim1$lane <- "stim1"
cells_stim1$type <- cells_stim1$marker

cells_stim2 <- Matrix::rowSums(cells_stim2@assays[["ADT"]]) %>% as.data.frame() 
cells_stim2$marker<- rownames(cells_stim2)
names(cells_stim2)[1] <- 'Counts'
cells_stim2$sample<- "cells_stim2"
cells_stim2$lane <- "stim2"
cells_stim2$type <- cells_stim2$marker

cells_stim3 <- Matrix::rowSums(cells_stim3@assays[["ADT"]]) %>% as.data.frame() 
cells_stim3$marker<- rownames(cells_stim3)
names(cells_stim3)[1] <- 'Counts'
cells_stim3$sample<- "cells_stim3"
cells_stim3$lane <- "stim3"
cells_stim3$type <- cells_stim3$marker

cells_stim4 <- Matrix::rowSums(cells_stim4@assays[["ADT"]]) %>% as.data.frame() 
cells_stim4$marker<- rownames(cells_stim4)
names(cells_stim4)[1] <- 'Counts'
cells_stim4$sample<- "cells_stim4"
cells_stim4$lane <- "stim4"
cells_stim4$type <- cells_stim4$marker

cells_stim5 <- Matrix::rowSums(cells_stim5@assays[["ADT"]]) %>% as.data.frame() 
cells_stim5$marker<- rownames(cells_stim5)
names(cells_stim5)[1] <- 'Counts'
cells_stim5$sample<- "cells_stim5"
cells_stim5$lane <- "stim5"
cells_stim5$type <- cells_stim5$marker

cells_unstim <- Matrix::rowSums(cells_unstim@assays[["ADT"]]) %>% as.data.frame() 
cells_unstim$marker<- rownames(cells_unstim)
names(cells_unstim)[1] <- 'Counts'
cells_unstim$sample<- "cells_unstim"
cells_unstim$lane <- "unstim"
cells_unstim$type <- cells_unstim$marker

all_stim1 <- rbind(empty_droplets_stim1,  cells_stim1)
all_stim2 <- rbind(empty_droplets_stim2, cells_stim2)
all_stim3 <- rbind(empty_droplets_stim3, cells_stim3)
all_stim4 <- rbind(empty_droplets_stim4, cells_stim4)
all_stim5 <- rbind(empty_droplets_stim5, cells_stim5)
all_unstim <- rbind(empty_droplets_unstim, cells_unstim)

Marker_order <- empty_droplets_stim1  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim1 <- Marker_order$marker

Marker_order <- empty_droplets_stim2  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim2 <- Marker_order$marker

Marker_order <- empty_droplets_stim3  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim3 <- Marker_order$marker

Marker_order <- empty_droplets_stim4  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim4 <- Marker_order$marker

Marker_order <- empty_droplets_stim5  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_stim5 <- Marker_order$marker

Marker_order <- empty_droplets_unstim  %>% group_by(marker) %>% 
  summarise(Frequency = sum(Counts))
Marker_order <- Marker_order[order(Marker_order$Frequency),]
order_unstim <- Marker_order$marker

all_stim1$marker <- factor(all_stim1$marker, levels = order_stim1)
all_stim2$marker <- factor(all_stim2$marker, levels = order_stim2)
all_stim3$marker <- factor(all_stim3$marker, levels = order_stim3)
all_stim4$marker <- factor(all_stim4$marker, levels = order_stim4)
all_stim5$marker <- factor(all_stim5$marker, levels = order_stim5)
all_unstim$marker <- factor(all_unstim$marker, levels = order_unstim)

Sti1 <- ggplot(all_stim1, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") 
Sti2 <- ggplot(all_stim2, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") 
Sti3 <- ggplot(all_stim3, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT")
Sti4 <- ggplot(all_stim4, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") 
Sti5 <- ggplot(all_stim5, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT") 
Uns <- ggplot(all_unstim, aes( x = "UMI", y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_light() + labs( x = "UMI", y="UMI Count (x 10^+6)", title = "ADT")
Sti1 +Sti2+ Sti3 +Sti4+ Sti5 + Uns

Sti1 <- ggplot(all_stim1, aes( x = marker, y = Counts/10^6, fill = sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("darkorchid",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim1") + theme(legend.position="top")
Sti2 <- ggplot(all_stim2, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("green",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim2") + theme(legend.position="top")
Sti3 <- ggplot(all_stim3, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("orange",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim3") + theme(legend.position="top")
Sti4 <- ggplot(all_stim4, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("orange",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim4") + theme(legend.position="top")
Sti5 <- ggplot(all_stim5, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("orange",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Stim5") + theme(legend.position="top")
Uns <- ggplot(all_unstim, aes( x = marker, y = Counts/10^6, fill=sample)) + geom_bar(stat = "identity", position = "stack", na.rm = TRUE)  + 
  theme_bw() +  scale_fill_manual(values = c("blue",  "gray48")) + 
  labs( x = "", y="UMI Count (x 10^+6)", title = "Unstim") + theme(legend.position="top")
Sti1 + coord_flip() | Sti2 + coord_flip() | Sti3 + coord_flip() |Sti4 + coord_flip()|Sti5 + coord_flip()| Uns + coord_flip()

#-----------------------------------------objects moving forward-----------------------------------------
cells_stim1 <- subset(x = stim1, subset = rank <= 50234)
cells_stim2 <- subset(x = stim2, subset = rank <= 43578)
cells_stim3 <- subset(x = stim3, subset = rank <= 27211) #based off of number of cells listed on cell ranger HTML summary 
cells_stim4 <- subset(x = stim4, subset = rank <= 41298) 
cells_stim5 <- subset(x = stim5, subset = rank <= 35376) 
cells_unstim <- subset(x = unstim, subset = rank <= 25356)

#-----------------------------------------demultiplex using HTODemux--------------------------------------------------
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cells_stim1 <- NormalizeData(cells_stim1, assay = "HTO", normalization.method = "CLR")
cells_stim2 <- NormalizeData(cells_stim2, assay = "HTO", normalization.method = "CLR")
cells_stim3 <- NormalizeData(cells_stim3, assay = "HTO", normalization.method = "CLR")
cells_stim4 <- NormalizeData(cells_stim4, assay = "HTO", normalization.method = "CLR")
cells_stim5 <- NormalizeData(cells_stim5, assay = "HTO", normalization.method = "CLR")
cells_unstim <- NormalizeData(cells_unstim, assay = "HTO", normalization.method = "CLR")
#We perform a k-medoid clustering on the normalized HTO values, which initially separates cells into K(# of samples)+1 clusters.
#We calculate a ‘negative’ distribution for HTO. For each HTO, we use the cluster with the lowest average value as the negative group.
#For each HTO, we fit a negative binomial distribution to the negative cluster. We use the 0.99 quantile of this distribution as a threshold.
#Based on these thresholds, each cell is classified as positive or negative for each HTO.
#Cells that are positive for more than one HTOs are annotated as doublets.

#positive.quantile - The quantile of inferred 'negative' distribution for each hashtag - over which the cell is considered 'positive'. Default is 0.99
#Determine how less stringent demultiplex would look by the highest HTO signal at or above a set UMI threshold
cells_stim1 <- HTODemux(cells_stim1, assay = "HTO", positive.quantile = 0.99)
cells_stim2 <- HTODemux(cells_stim2, assay = "HTO", positive.quantile = 0.99)
cells_stim3 <- HTODemux(cells_stim3, assay = "HTO", positive.quantile = 0.99)
cells_stim4 <- HTODemux(cells_stim4, assay = "HTO", positive.quantile = 0.99)
cells_stim5 <- HTODemux(cells_stim5, assay = "HTO", positive.quantile = 0.99)
cells_unstim <- HTODemux(cells_unstim, assay = "HTO", positive.quantile = 0.99)

HTO_levels <- c("Negative", "Subject1-7moPostVax","Subject1-1moPostBoost","Subject1-PostBreakThru","Subject4-PostBreakThru",
                "Subject4-7moPostVax","Subject4-1moPostBoost","Subject10-7moPostVax","Subject10-1moPostBoost",
                "Subject11-7moPostVax","Subject11-1moPostBoost","Subject5-PostBreakThru", "Doublet")

##Stim1
table(cells_stim1$HTO_classification.global) 
table(cells_stim1$HTO_maxID) 
cells_stim1$hash.ID <- factor(cells_stim1$hash.ID, levels = HTO_levels)
Idents(cells_stim1) <- "HTO_maxID"
s.1 <- VlnPlot(cells_stim1, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "stim") +  NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
s.1.1 <- RidgePlot(cells_stim1, assay = "HTO", features = rownames(cells_stim1[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 4)

##stim.2
table(cells_stim2$HTO_classification.global)
table(cells_stim2$HTO_maxID)
cells_stim2$hash.ID <- factor(cells_stim2$hash.ID, levels = HTO_levels)
Idents(cells_stim2) <- "HTO_maxID"
s.2<- VlnPlot(cells_stim2, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "tetramer") + NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
s.2.1 <- RidgePlot(cells_stim2, assay = "HTO", features = rownames(cells_stim2[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 1)

##stim.3
table(cells_stim3$HTO_classification.global)
table(cells_stim3$HTO_maxID)
cells_stim3$hash.ID <- factor(cells_stim3$hash.ID, levels = HTO_levels)
Idents(cells_stim3) <- "HTO_maxID"
s.3<- VlnPlot(cells_stim3, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "tetramer") + NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
s.3.1 <- RidgePlot(cells_stim3, assay = "HTO", features = rownames(cells_stim3[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 1)

##stim.4
table(cells_stim4$HTO_classification.global)
table(cells_stim4$HTO_maxID)
cells_stim4$hash.ID <- factor(cells_stim4$hash.ID, levels = HTO_levels)
Idents(cells_stim4) <- "HTO_maxID"
s.4<- VlnPlot(cells_stim4, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "tetramer") + NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
s.4.1 <- RidgePlot(cells_stim4, assay = "HTO", features = rownames(cells_stim4[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 1)

##stim.5
table(cells_stim5$HTO_classification.global)
table(cells_stim5$hash.ID)
cells_stim5$hash.ID <- factor(cells_stim5$hash.ID, levels = HTO_levels)
Idents(cells_stim5) <- "HTO_maxID"
s.5<- VlnPlot(cells_stim5, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "tetramer") + NoLegend() + 
  geom_hline(yintercept = 10)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
s.5.1 <- RidgePlot(cells_stim5, assay = "HTO", features = rownames(cells_stim5[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 4)

##unstim
table(cells_unstim$HTO_classification.global)
table(cells_unstim$HTO_maxID)

Idents(cells_unstim) <- "HTO_maxID"
u<- VlnPlot(cells_unstim, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "unstim") + NoLegend() + 
  geom_hline(yintercept = 10)+ scale_y_continuous(limits = c(1,100000), trans = "log10")
u.1 <- RidgePlot(cells_unstim, assay = "HTO", features = rownames(cells_unstim[["HTO"]]), group.by="hash.ID", combine=TRUE, ncol = 1)

#FEATURE SCATTER PLOT OF HTOS 
s.1 | s.2 | s.3 | s.4 | s.5 | u.1 #violin plots w nCount_HTO
s.1.1 | s.2.1 | s.3.1 | s.4.1 | s.5.1 | u.1.1 #ridge plots

#-----------------------------------------demultiplex using less stringent HTODemux--------------------------------------------------
## Determine how less stringent demultiplex would look with HTO_maxID at or above a set UMI threshold of ??????
cells_stim1 <- subset(cells_stim1, subset = nCount_HTO >= 10)
cells_stim2 <- subset(cells_stim2, subset = nCount_HTO >= 10)
cells_stim3 <- subset(cells_stim3, subset= nCount_HTO >= 10)
cells_stim4 <- subset(cells_stim4, subset= nCount_HTO >= 10)
cells_stim5 <- subset(cells_stim5, subset= nCount_HTO >= 10)
cells_unstim <- subset(cells_unstim, subset = nCount_HTO >= 10)

stim_ID1 <- cells_stim1$HTO_classification.global
stim_ID1[cells_stim1@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim1$stim_ID1 <- factor(stim_ID1)

stim_ID2 <- cells_stim2$HTO_classification.global
stim_ID2[cells_stim2@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim2$stim_ID2 <- factor(stim_ID2)

stim_ID3 <- cells_stim3$HTO_classification.global
stim_ID3[cells_stim3@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim3$stim_ID3 <- factor(stim_ID3)

stim_ID4 <- cells_stim4$HTO_classification.global
stim_ID4[cells_stim4@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim4$stim_ID4 <- factor(stim_ID4)

stim_ID5 <- cells_stim5$HTO_classification.global
stim_ID5[cells_stim5@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_stim5$stim_ID5 <- factor(stim_ID5)

unstim_ID <- cells_unstim$HTO_classification.global
unstim_ID[cells_unstim@meta.data$HTO_classification.global == "Negative"] <- "Singlet"
cells_unstim$unstim_ID <- factor(unstim_ID)

#look at numbers of singlets and doublets
s.1 <- cells_stim1@meta.data %>% 
  ggplot(aes(x= stim_ID1)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim1 NCells") + ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.2 <- cells_stim2@meta.data %>% 
  ggplot(aes(x=stim_ID2)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim2 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.3 <- cells_stim3@meta.data %>% 
  ggplot(aes(x=stim_ID3)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim3 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.4 <- cells_stim4@meta.data %>% 
  ggplot(aes(x=stim_ID4)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim4 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.5 <- cells_stim5@meta.data %>% 
  ggplot(aes(x=stim_ID5)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim5 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

u.1 <-cells_unstim@meta.data %>% 
  ggplot(aes(x=unstim_ID)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.1 | s.2 
s.3 | s.4  
s.5 | u.1
#-----------------------------------------REMOVE UNUSED OBJECTS-----------------------------------------
remove(stim1.out,stim2.out,stim3.out,stim4.out, stim5.out, unstim.out, adt_stim1,adt_stim2,adt_stim3,adt_stim4,adt_stim5,adt_unstim,
       stim1, stim2, stim3,stim4,stim5, unstim, 
       adt_assay_stim1, adt_assay_stim2, adt_assay_stim3,adt_assay_stim4,adt_assay_stim5, adt_assay_unstim,
       HTO_assay_stim1, HTO_assay_stim2, HTO_assay_stim3,HTO_assay_stim4,HTO_assay_stim5, HTO_assay_unstim)

#-----------------------------------------stim1-----------------------------------------
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim1@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim1@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim1@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim1 <- AddMetaData(object = cells_stim1, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim1, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "stim")
FeatureScatter(object = cells_stim1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "seurat_clusters") +  geom_hline(yintercept = c(400, 7000)) + labs(title = "stim")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_stim1) <- "HTO_classification.global"
RidgePlot(cells_stim1, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim1 <- subset(x = cells_stim1, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim1 <- NormalizeData(cells_stim1, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim1 <- FindVariableFeatures(cells_stim1, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim1))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim1 <- ScaleData(cells_stim1, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim1 <- RunPCA(cells_stim1, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim1, reduction = "pca")
DimHeatmap(object = cells_stim1, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)
cells_stim1 <- FindNeighbors(cells_stim1, dims = 1:30)
cells_stim1 <- FindClusters(cells_stim1, resolution = 0.8, verbose = FALSE)
cells_stim1 <- RunUMAP(cells_stim1, dims = 1:30)
DimPlot(cells_stim1, label = TRUE, group.by = "seurat_clusters")
s1 <- DimPlot(cells_stim1, group.by = "seurat_clusters",label = TRUE) + NoLegend()

# Normalize ADT data,
cells_stim1 <- NormalizeData(cells_stim1, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim1 <- ScaleData(cells_stim1, assay = "ADT")

#-----------------------------------------stim2-----------------------------------------
DefaultAssay(cells_stim2) <- "RNA"
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim2@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim2@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim2@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim2 <- AddMetaData(object = cells_stim2, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident") 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim2, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "stim2")
FeatureScatter(object = cells_stim2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification.global") +  geom_hline(yintercept = c(400, 7000)) + labs(title = "stim2")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_stim2) <- "HTO_classification.global"
RidgePlot(cells_stim2, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim2 <- subset(x = cells_stim2, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim2 <- NormalizeData(cells_stim2, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim2 <- FindVariableFeatures(cells_stim2, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim2))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim2 <- ScaleData(cells_stim2, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim2 <- RunPCA(cells_stim2, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim2, reduction = "pca")
DimHeatmap(object = cells_stim2, reduction = "pca", cells = 500, balanced = TRUE)

cells_stim2 <- FindNeighbors(cells_stim2, dims = 1:30)
cells_stim2 <- FindClusters(cells_stim2, resolution = 0.8, verbose = FALSE)
cells_stim2 <- RunUMAP(cells_stim2, dims = 1:30)
s2<- DimPlot(cells_stim2, label = TRUE) + NoLegend()

# Normalize ADT data,
cells_stim2 <- NormalizeData(cells_stim2, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim2 <- ScaleData(cells_stim2, assay = "ADT")
#-----------------------------------------stim3-----------------------------------------
DefaultAssay(cells_stim3) <- "RNA"
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim3@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim3@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim3@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim3 <- AddMetaData(object = cells_stim3, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident") 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim3, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "whole PBMC")
FeatureScatter(object = cells_stim3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +  geom_hline(yintercept = c(400, 7000)) + labs(title = "stim3")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_stim3) <- "HTO_classification.global"
RidgePlot(cells_stim3, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim3 <- subset(x = cells_stim3, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim3 <- NormalizeData(cells_stim3, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim3 <- FindVariableFeatures(cells_stim3, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim3))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim3 <- ScaleData(cells_stim3, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim3 <- RunPCA(cells_stim3, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim3, reduction = "pca")
DimHeatmap(object = cells_stim3, reduction = "pca", cells = 500, balanced = TRUE)

cells_stim3 <- FindNeighbors(cells_stim3, dims = 1:30)
cells_stim3 <- FindClusters(cells_stim3, resolution = 0.8, verbose = FALSE)
cells_stim3 <- RunUMAP(cells_stim3, dims = 1:30)
s3 <- DimPlot(cells_stim3, label = TRUE)  + NoLegend()

# Normalize ADT data,
cells_stim3 <- NormalizeData(cells_stim3, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim3 <- ScaleData(cells_stim3, assay = "ADT")
#-----------------------------------------stim4-----------------------------------------
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim4@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim4@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim4@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim4 <- AddMetaData(object = cells_stim4, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim4, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim4, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "stim")
FeatureScatter(object = cells_stim4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification.global") +  geom_hline(yintercept = c(400, 7000)) + labs(title = "stim")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_stim4) <- "HTO_classification.global"
RidgePlot(cells_stim4, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim4 <- subset(x = cells_stim4, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim4 <- NormalizeData(cells_stim4, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim4 <- FindVariableFeatures(cells_stim4, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim4))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim4 <- ScaleData(cells_stim4, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim4 <- RunPCA(cells_stim4, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim4, reduction = "pca")
DimHeatmap(object = cells_stim4, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)
cells_stim4 <- FindNeighbors(cells_stim4, dims = 1:30)
cells_stim4 <- FindClusters(cells_stim4, resolution = 0.8, verbose = FALSE)
cells_stim4 <- RunUMAP(cells_stim4, dims = 1:30)
DimPlot(cells_stim4, label = TRUE, group.by = "seurat_clusters")
s4 <- DimPlot(cells_stim4, group.by = "seurat_clusters",label = TRUE) + NoLegend()

# Normalize ADT data,
cells_stim4 <- NormalizeData(cells_stim4, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim4 <- ScaleData(cells_stim4, assay = "ADT")
#-----------------------------------------stim5-----------------------------------------
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_stim5@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_stim5@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_stim5@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_stim5 <- AddMetaData(object = cells_stim5, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_stim5, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_stim5, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "stim")
FeatureScatter(object = cells_stim5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "HTO_classification.global") +  geom_hline(yintercept = c(400, 7000)) + labs(title = "stim")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_stim5) <- "HTO_classification.global"
RidgePlot(cells_stim5, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_stim5 <- subset(x = cells_stim5, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

#perform visualization and clustering steps 
#Normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
cells_stim5 <- NormalizeData(cells_stim5, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins
#to calculate z-score for dispersion within each bin
cells_stim5 <- FindVariableFeatures(cells_stim5, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
head(x = HVFInfo(object = cells_stim5))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
cells_stim5 <- ScaleData(cells_stim5, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
cells_stim5 <- RunPCA(cells_stim5, verbose = FALSE)

# Examine and visualize PCA results a few different ways
DimPlot(object = cells_stim5, reduction = "pca")
DimHeatmap(object = cells_stim5, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)
cells_stim5 <- FindNeighbors(cells_stim5, dims = 1:30)
cells_stim5 <- FindClusters(cells_stim5, resolution = 0.8, verbose = FALSE)
cells_stim5 <- RunUMAP(cells_stim5, dims = 1:30)
DimPlot(cells_stim5, label = TRUE, group.by = "seurat_clusters")
s5 <- DimPlot(cells_stim5, group.by = "seurat_clusters",label = TRUE) + NoLegend()

# Normalize ADT data,
cells_stim5 <- NormalizeData(cells_stim5, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_stim5 <- ScaleData(cells_stim5, assay = "ADT")
#--------------DIM plots across three strim lanes---------------------------------
s1 | s2 | s3 | s4 | s5

#-----------------------------------------unstim-----------------------------------------
DefaultAssay(cells_unstim) <- "RNA"
#The % of UMI mapping to MT-genes is a common scRNA-seq QC metric 
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(cells_unstim@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(cells_unstim@assays[["RNA"]][mito.genes, ])/Matrix::colSums(cells_unstim@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
cells_unstim <- AddMetaData(object = cells_unstim, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = cells_unstim, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident") 

#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = cells_unstim, feature1 = "nCount_RNA", feature2 = "percent.mito",group.by = "orig.ident" ) + geom_hline(yintercept = 0.08) + labs(title = "unstim")
FeatureScatter(object = cells_unstim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +  geom_hline(yintercept = c(200, 7000)) + labs(title = "unstim")

#Do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(cells_unstim) <- "hash.ID"
RidgePlot(cells_unstim, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

# We filter out cells that have unique gene counts (nFeature_RNA) over 4,000 or less than 200 Note that > and < are used to define a 'gate' 
# -Inf and Inf should be used if you don't want a lower or upper threshold.
cells_unstim <- subset(x = cells_unstim, subset = nFeature_RNA > 400 & nFeature_RNA < 7000 & percent.mito >  -Inf & percent.mito < 0.08 ) 

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
DimPlot(cells_unstim, label = FALSE)

# Normalize ADT data,
cells_unstim <- NormalizeData(cells_unstim, normalization.method = "CLR", margin = 2, assay = "ADT") 
cells_unstim <- ScaleData(cells_unstim, assay = "ADT")

#-----------------------------saving objects---------------------------
setwd("/Volumes/heratr01lab/heratr01labspace/ECCITE-seq/Differential_Spike_Priming_Studies/Experiments/20220927_ECCITEseq_AIM")
saveRDS(cells_stim1, file = "20221021_stim1.rds")
saveRDS(cells_stim2, file = "20221021_stim2.rds")
saveRDS(cells_stim3, file = "20221021_stim3.rds")
saveRDS(cells_stim4, file = "20221021_stim4.rds")
saveRDS(cells_stim5, file = "20221021_stim5.rds")
saveRDS(cells_unstim, file = "20221021_unstim.rds")

cells_stim1 <- readRDS(file = "20221021_stim1.rds")
cells_stim2<- readRDS( file = "20221021_stim2.rds")
cells_stim3<- readRDS(file = "20221021_stim3.rds")
cells_stim4<- readRDS(file = "20221021_stim4.rds")
cells_stim5<- readRDS( file = "20221021_stim5.rds")
cells_unstim<- readRDS(file = "20221021_unstim.rds")
#-----------------------------------------filter out doublets--------------------------------------------------
library(scDblFinder) #new doublet rate since using the HT kit...assumed to be 0.4% per thousand cells 
DefaultAssay(cells_stim1) <- "RNA"
DefaultAssay(cells_stim2) <- "RNA"
DefaultAssay(cells_stim3) <- "RNA"
DefaultAssay(cells_stim4) <- "RNA"
DefaultAssay(cells_stim5) <- "RNA"
DefaultAssay(cells_unstim) <- "RNA"

sce.stim1 <- as.SingleCellExperiment(cells_stim1)
sce.stim1 <- scDblFinder(sce.stim1, dbr = (0.004 * ncol(sce.stim1)/1000), samples="HTO_maxID")
table(sce.stim1$scDblFinder.class, sce.stim1$HTO_maxID)

sce.stim2 <- as.SingleCellExperiment(cells_stim2)
sce.stim2 <- scDblFinder(sce.stim2, dbr = (0.004 * ncol(sce.stim2)/1000), samples="HTO_maxID")
table(sce.stim2$scDblFinder.class, sce.stim2$HTO_maxID)

sce.stim3 <- as.SingleCellExperiment(cells_stim3)
sce.stim3 <- scDblFinder(sce.stim3, dbr = (0.004 * ncol(sce.stim3)/1000), samples="HTO_maxID")
table(sce.stim3$scDblFinder.class, sce.stim3$HTO_maxID)

sce.stim4 <- as.SingleCellExperiment(cells_stim4)
sce.stim4 <- scDblFinder(sce.stim4, dbr = (0.004 * ncol(sce.stim4)/1000), samples="HTO_maxID")
table(sce.stim4$scDblFinder.class, sce.stim4$HTO_maxID)

sce.stim5 <- as.SingleCellExperiment(cells_stim5)
sce.stim5 <- scDblFinder(sce.stim5, dbr = (0.004 * ncol(sce.stim5)/1000), samples="HTO_maxID")
table(sce.stim5$scDblFinder.class, sce.stim5$HTO_maxID)

sce.uns <- as.SingleCellExperiment(cells_unstim)
sce.uns <- scDblFinder(sce.uns, dbr = (0.004 * ncol(sce.uns)/1000), samples="HTO_maxID")
table(sce.uns$scDblFinder.class, sce.uns$HTO_maxID)

identical(colnames(cells_stim1),colnames(sce.stim1))
identical(colnames(cells_stim2),colnames(sce.stim2))
identical(colnames(cells_stim3),colnames(sce.stim3))
identical(colnames(cells_stim4),colnames(sce.stim4))
identical(colnames(cells_stim5),colnames(sce.stim5))
identical(colnames(cells_unstim),colnames(sce.uns))

cells_stim1$scDblFinder.class <- sce.stim1$scDblFinder.class
cells_stim1$scDblFinder.score <- sce.stim1$scDblFinder.score
cells_stim1$scDblFinder.weighted <- sce.stim1$scDblFinder.weighted
rm(sce.stim1)

cells_stim2$scDblFinder.class <- sce.stim2$scDblFinder.class
cells_stim2$scDblFinder.score <- sce.stim2$scDblFinder.score
cells_stim2$scDblFinder.weighted <- sce.stim2$scDblFinder.weighted
rm(sce.stim2)

cells_stim3$scDblFinder.class <- sce.stim3$scDblFinder.class
cells_stim3$scDblFinder.score <- sce.stim3$scDblFinder.score
cells_stim3$scDblFinder.weighted <- sce.stim3$scDblFinder.weighted
rm(sce.stim3)

cells_stim4$scDblFinder.class <- sce.stim4$scDblFinder.class
cells_stim4$scDblFinder.score <- sce.stim4$scDblFinder.score
cells_stim4$scDblFinder.weighted <- sce.stim4$scDblFinder.weighted
rm(sce.stim4)

cells_stim5$scDblFinder.class <- sce.stim5$scDblFinder.class
cells_stim5$scDblFinder.score <- sce.stim5$scDblFinder.score
cells_stim5$scDblFinder.weighted <- sce.stim5$scDblFinder.weighted
rm(sce.stim5)

cells_unstim$scDblFinder.class <- sce.uns$scDblFinder.class
cells_unstim$scDblFinder.score <- sce.uns$scDblFinder.score
cells_unstim$scDblFinder.weighted <- sce.uns$scDblFinder.weighted
rm(sce.uns)

cells_stim1$scDblFinder.class <- factor(cells_stim1$scDblFinder.class, levels=c("doublet", "singlet"))
cells_stim2$scDblFinder.class <- factor(cells_stim2$scDblFinder.class, levels=c("doublet", "singlet"))
cells_stim3$scDblFinder.class <- factor(cells_stim3$scDblFinder.class, levels=c("doublet", "singlet"))
cells_stim4$scDblFinder.class <- factor(cells_stim4$scDblFinder.class, levels=c("doublet", "singlet"))
cells_stim4$scDblFinder.class <- factor(cells_stim4$scDblFinder.class, levels=c("doublet", "singlet"))
cells_unstim$scDblFinder.class <- factor(cells_unstim$scDblFinder.class, levels=c("doublet", "singlet"))

#look at numbers of singlets and doublets using scDblFinder
s.1 <- cells_stim1@meta.data %>% 
  ggplot(aes(x= scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim1 NCells") + ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.2 <- cells_stim2@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim2 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.3 <-cells_stim3@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim3 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.4 <-cells_stim4@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim4 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.5 <-cells_stim5@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim5 NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

u.1 <-cells_unstim@meta.data %>% 
  ggplot(aes(x=scDblFinder.class)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells")+ ylim(c(0,40000)) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

s.1 | s.2 
s.3 | s.4  
s.5 | u.1


#-----------------------------------Filtering out doublets?-------------------------------------------
table(cells_stim1$scDblFinder.class)
table(cells_stim2$scDblFinder.class)
table(cells_stim3$scDblFinder.class)
table(cells_stim4$scDblFinder.class)
table(cells_stim5$scDblFinder.class)
table(cells_unstim$scDblFinder.class)

#CD19 vs CD3 gene expression before filtering out doublets
s.1 <- FeatureScatter(cells_stim1, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() + labs(title="stim1") #+ geom_vline(xintercept = 0.75) +  geom_hline(yintercept = 0.75)
s.2 <- FeatureScatter(cells_stim2, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() + labs(title="stim2") #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)
s.3 <- FeatureScatter(cells_stim3, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() + labs(title="stim3") #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)
s.4 <- FeatureScatter(cells_stim4, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() + labs(title="stim4") #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)
s.5 <- FeatureScatter(cells_stim5, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() + labs(title="stim5") #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)
u.1 <- FeatureScatter(cells_unstim, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + labs(title="unstim")# + NoLegend() #+geom_vline(xintercept = 0.75)+ geom_hline(yintercept = 0.75)

s.1 | s.2 |s.3 | s.4  |s.5 | u.1

#subset to on singlets 
stim1_singlet <- subset(cells_stim1, subset=scDblFinder.class=="singlet")
stim2_singlet <- subset(cells_stim2, subset=scDblFinder.class=="singlet")
stim3_singlet <- subset(cells_stim3, subset=scDblFinder.class=="singlet")
stim4_singlet <- subset(cells_stim4, subset=scDblFinder.class=="singlet")
stim5_singlet <- subset(cells_stim5, subset=scDblFinder.class=="singlet")
uns_singlet <- subset(cells_unstim, subset=scDblFinder.class=="singlet")

#subset singlets round 2 by excluding cells that are double positive for CD19 and CD3E genes 
s.1 <- FeatureScatter(cells_stim1, feature1 = "CD3E", feature2 = "CD19", group.by = "stim_ID1") + NoLegend() + labs(title="stim1") + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim1") + NoLegend()
s.2 <- FeatureScatter(cells_stim2, feature1 = "CD3E", feature2 = "CD19", group.by = "stim_ID2") + NoLegend() + labs(title="stim2") + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim2") + NoLegend()
s.3 <- FeatureScatter(cells_stim3, feature1 = "CD3E", feature2 = "CD19", group.by = "stim_ID3") + NoLegend() + labs(title="stim3")+ geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim3") + NoLegend()
s.4 <- FeatureScatter(cells_stim4, feature1 = "CD3E", feature2 = "CD19", group.by = "stim_ID4") + NoLegend() + labs(title="stim4")+ geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim4") + NoLegend()
s.5 <- FeatureScatter(cells_stim5, feature1 = "CD3E", feature2 = "CD19", group.by = "stim_ID5") + NoLegend() + labs(title="stim5") + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim5") + NoLegend()
u.1 <- FeatureScatter(cells_unstim, feature1 = "CD3E", feature2 = "CD19", group.by = "unstim_ID") + labs(title="unstim")+geom_vline(xintercept = 0.1)+ geom_hline(yintercept = 0.1) + labs(title = "unstim")

s.1 | s.2 |s.3 | s.4  |s.5 | u.1

stim1_singlet.2 <- subset(stim1_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)
stim2_singlet.2 <- subset(stim2_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)
stim3_singlet.2 <- subset(stim3_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)
stim4_singlet.2 <- subset(stim4_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)
stim5_singlet.2 <- subset(stim5_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)
uns_singlet.2 <- subset(uns_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1| CD3E < 0.1 & CD19 < 0.1 )

FeatureScatter(stim1_singlet.2, feature1 = "CD3E", feature2 = "CD19", group.by = "HTO_maxID") + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1) + labs(title = "stim")

#look at number of singlets
s.1 <- stim1_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim1 NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

s.2 <- stim2_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim2 NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

s.3 <- stim3_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim3 NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

s.4 <- stim4_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim4 NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

s.5 <- stim5_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Stim5 NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

u.1 <-uns_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle(" Unstim NCells")+
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

s.1 | s.2 |s.3 | s.4  |s.5 | u.1

#------------------------------------SAVING final objects--------------------------------------------------
stim1_singlet.2@meta.data <- stim1_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
stim1_singlet.2@meta.data <- stim1_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
stim1_singlet.2@meta.data <- stim1_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no",
    endsWith(subject, "Subject4") ~ "no",
    endsWith(subject, "Subject5") ~ "no"
  ))
stim1_singlet.2$timepoint <- factor(stim1_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

stim2_singlet.2@meta.data <- stim2_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
stim2_singlet.2@meta.data <- stim2_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
stim2_singlet.2@meta.data <- stim2_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no",
    endsWith(subject, "Subject4") ~ "no",
    endsWith(subject, "Subject5") ~ "no"
  ))
stim2_singlet.2$timepoint <- factor(stim2_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

stim3_singlet.2@meta.data <- stim3_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
stim3_singlet.2@meta.data <- stim3_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
stim3_singlet.2@meta.data <- stim3_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no",
    endsWith(subject, "Subject4") ~ "no",
    endsWith(subject, "Subject5") ~ "no"
  ))
stim3_singlet.2$timepoint <- factor(stim3_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

stim4_singlet.2@meta.data <- stim4_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
stim4_singlet.2@meta.data <- stim4_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
stim4_singlet.2@meta.data <- stim4_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no",
    endsWith(subject, "Subject4") ~ "no",
    endsWith(subject, "Subject5") ~ "no"
  ))
stim4_singlet.2$timepoint <- factor(stim4_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

stim5_singlet.2@meta.data <- stim5_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
stim5_singlet.2@meta.data <- stim5_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
stim5_singlet.2@meta.data <- stim5_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no",
    endsWith(subject, "Subject4") ~ "no",
    endsWith(subject, "Subject5") ~ "no"
  ))
stim5_singlet.2$timepoint <- factor(stim5_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

uns_singlet.2@meta.data <- uns_singlet.2@meta.data %>% separate(HTO_maxID, c("subject", "subject1", "timepoint"), sep = "-", remove = FALSE)
uns_singlet.2@meta.data <- uns_singlet.2@meta.data %>% unite("subject", subject:subject1, remove=TRUE, sep = "-")
uns_singlet.2@meta.data <- uns_singlet.2@meta.data  %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject11")  ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",
    endsWith(subject, "Subject1")  ~ "no"))
uns_singlet.2$timepoint <- factor(uns_singlet.2$timepoint, levels=c("7moPostVax","1moPostBoost", "PostBreakThru"))

setwd("/Volumes/heratr01lab/heratr01labspace/ECCITE-seq/Differential_Spike_Priming_Studies/Experiments/20220927_ECCITEseq_AIM")
saveRDS(stim1_singlet.2, file = "20221021_stim1.rds")
saveRDS(stim2_singlet.2, file = "20221021_stim2.rds")
saveRDS(stim3_singlet.2, file = "20221021_stim3.rds")
saveRDS(stim4_singlet.2, file = "20221021_stim4.rds")
saveRDS(stim5_singlet.2, file = "20221021_stim5.rds")
saveRDS(uns_singlet.2, file = "20221021_unstim.rds")
