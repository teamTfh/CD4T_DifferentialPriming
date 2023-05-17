library(pheatmap)
library(SingleCellExperiment)
library(ggplot2)
library(gplots)
library(tidyverse)
library(tidyr)
library(dplyr)
#devtools::install_github("ncborcherding/scRepertoire@dev", force = TRUE)
#BiocManager::install("scRepertoire")
library(scRepertoire)
library(readxl)
library(DESeq2)
library(SummarizedExperiment)
library(GSEABase)
library(GSVA)
library(ggbiplot)
#BiocManager::install("glmGamPoi")
library(glmGamPoi)
library(Seurat)
#BiocManager::install("scuttle")
library(scuttle)
#BiocManager::install("limma")
library(limma)

sessionInfo()
setwd("~/S3data/Saved_Objects")
#-------colors------
pre_vax <- "#fdc472"
post_vax <- "#e59628"
pre_infect <- "#b5b2f1"
post_infect <- "#6667a9"
vax <- "#ffc471"
infect <- "#b5b2f1"
breakthru <- "#D96F6F"
colorblind_vector <- colorRampPalette((c("#0D0887FF", "#0D0887FF", 
                                         "#0D0887FF", "#47039FFF", "#BD3786FF", "#BD3786FF",
                                         "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
                                                                  
#-------------------------Reading in pre-processed objects-------------- 
#20211122_exCITEseq_processing
stim_NOV <- readRDS("20211114_stim.rds")
unstim_NOV <- readRDS("20211114_unstim.rds")
stim_NOV$HTO_maxID <-  gsub('7moPostVax', 'Subject2_7moPostVax', stim_NOV$HTO_maxID)
stim_NOV$HTO_maxID <-  gsub('1moPostBoost', 'Subject2_1moPostBoost', stim_NOV$HTO_maxID)
unstim_NOV$HTO_maxID <-  gsub('7moPostVax', 'Subject2_7moPostVax', unstim_NOV$HTO_maxID)
unstim_NOV$HTO_maxID <-  gsub('1moPostBoost', 'Subject2_1moPostBoost', unstim_NOV$HTO_maxID)

#20211216_exCITEseq_processing
stim_DEC <- readRDS("20211203_stim.rds")
unstim_DEC <- readRDS("20211203_unstim.rds")
#View(stim_DEC@meta.data)

#20220208_exCITEseq_processing 
stim1_JAN <- readRDS("20220131_stim1.rds")
stim2_JAN <- readRDS("20220131_stim2.rds")
unstim_JAN <- readRDS("20220131_unstim.rds")
#View(stim1_JAN@meta.data)

#20220401_exCITEseq_processing
stim1_MARCH <- readRDS("20220320_stim1.rds")
stim2_MARCH <- readRDS("20220320_stim2.rds")
stim3_MARCH <- readRDS("20220320_stim3.rds")
unstim_MARCH <- readRDS("20220320_unstim.rds")
#View(stim3_MARCH@meta.data)

#20221005_exCITEseq_processing
stim1_SEPT <- readRDS("20221021_stim1.rds")
stim2_SEPT <- readRDS("20221021_stim2.rds")
stim3_SEPT <- readRDS("20221021_stim3.rds")
stim4_SEPT <- readRDS("20221021_stim4.rds")
stim5_SEPT <- readRDS("20221021_stim5.rds")
unstim_SEPT <- readRDS("20221021_unstim.rds")

stim_NOV@meta.data$orig.ident <- "stim_NOV"
unstim_NOV@meta.data$orig.ident <- "unstim_NOV"
stim_DEC@meta.data$orig.ident <- "stim_DEC"
unstim_DEC@meta.data$orig.ident <- "unstim_DEC"
stim1_JAN@meta.data$orig.ident <- "stim1_JAN"
stim2_JAN@meta.data$orig.ident <- "stim2_JAN"
unstim_JAN@meta.data$orig.ident <- "unstim_JAN"
stim1_MARCH@meta.data$orig.ident <- "stim1_MARCH"
stim2_MARCH@meta.data$orig.ident <- "stim2_MARCH"
stim3_MARCH@meta.data$orig.ident <- "stim3_MARCH"
unstim_MARCH@meta.data$orig.ident <- "unstim_MARCH"
stim1_SEPT@meta.data$orig.ident <- "stim1_SEPT"
stim2_SEPT@meta.data$orig.ident <- "stim2_SEPT"
stim3_SEPT@meta.data$orig.ident <- "stim3_SEPT"
stim4_SEPT@meta.data$orig.ident <- "stim4_SEPT"
stim5_SEPT@meta.data$orig.ident <- "stim5_SEPT"
unstim_SEPT@meta.data$orig.ident <- "unstim_SEPT"

stim_NOV@meta.data$sample <- "stim"
stim_DEC@meta.data$sample <- "stim"
stim1_JAN@meta.data$sample <- "stim"
stim2_JAN@meta.data$sample <- "stim"
stim1_MARCH@meta.data$sample <- "stim"
stim2_MARCH@meta.data$sample <- "stim"
stim3_MARCH@meta.data$sample <- "stim"
stim1_SEPT@meta.data$sample <- "stim"
stim2_SEPT@meta.data$sample <- "stim"
stim3_SEPT@meta.data$sample <-"stim"
stim4_SEPT@meta.data$sample <-"stim"
stim5_SEPT@meta.data$sample <- "stim"
unstim_NOV@meta.data$sample <- "unstim"
unstim_DEC@meta.data$sample <- "unstim"
unstim_JAN@meta.data$sample <- "unstim"
unstim_MARCH@meta.data$sample <- "unstim"
unstim_SEPT@meta.data$sample <- "unstim"

stim_NOV@meta.data$run <- "November"
stim_DEC@meta.data$run <- "December"
stim1_JAN@meta.data$run <- "January"
stim2_JAN@meta.data$run <- "January"
unstim_NOV@meta.data$run <- "November"
unstim_DEC@meta.data$run <- "December"
unstim_JAN@meta.data$run <- "January"
stim1_MARCH@meta.data$run <- "March"
stim2_MARCH@meta.data$run <- "March"
stim3_MARCH@meta.data$run <- "March"
unstim_MARCH@meta.data$run <- "March"
stim1_SEPT@meta.data$run <- "September"
stim2_SEPT@meta.data$run <- "September"
stim3_SEPT@meta.data$run <- "September"
stim4_SEPT@meta.data$run <- "September"
stim5_SEPT@meta.data$run <- "September"
unstim_SEPT@meta.data$run <- "September"

#-------------------------Integrating all AIM stim and unstim lanes--------------
#SOURCE:
#https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#using SCTransform and RPCA to integrate datasets quickly...default k.anchors = 5 for FindIntegrationFeatures
AIM <- merge(stim_NOV, c(stim_DEC, stim1_JAN, stim2_JAN,stim1_MARCH, stim2_MARCH, stim3_MARCH,stim1_SEPT,stim2_SEPT,stim3_SEPT,stim4_SEPT,stim5_SEPT,
                         unstim_NOV, unstim_DEC, unstim_JAN, unstim_MARCH,unstim_SEPT)) 
View(AIM)
AIM$orig.ident <- as.factor(AIM$orig.ident)
levels(AIM$orig.ident)
DefaultAssay(AIM) <- "RNA"

# split the dataset into a list of the 11 seurat objects
AIM <- SplitObject(AIM, split.by = "orig.ident")
remove(stim_NOV,stim_DEC, stim1_JAN, stim2_JAN,stim1_MARCH, stim2_MARCH, stim3_MARCH,stim1_SEPT,stim2_SEPT,stim3_SEPT,stim4_SEPT,stim5_SEPT,
       unstim_NOV, unstim_DEC, unstim_JAN, unstim_MARCH,unstim_SEPT)

# normalize using SCTransform and identify variable features for each dataset independently
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
AIM <- lapply(X = AIM, FUN = SCTransform, method = "glmGamPoi")

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
#Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = AIM, nfeatures = 3000)
AIM <- PrepSCTIntegration(object.list = AIM, anchor.features = features)
AIM <- lapply(X = AIM, FUN = RunPCA, features = features)

DefaultAssay(AIM$stim_NOV) #should be SCT

#We then identify anchors using the FindIntegrationAnchors() function, which takes a list of 
#Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
immune.anchors <- FindIntegrationAnchors(object.list = AIM, normalization.method = "SCT", 
                                         reference = c(1,13),  #reference datasets are form november -- Stim and Unstim 
                                         anchor.features = features, dims = 1:30, 
                                         reduction = "rpca", k.anchor = 20)

saveRDS(features, file = "20221023_integration.features.rds")
saveRDS(immune.anchors, file = "20221023_immune.anchors.rds")
saveRDS(AIM, file = "20221023_AIM.SplitObject.rds")
rm(AIM)
AIM.integrated.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

AIM.integrated.sct <- RunPCA(AIM.integrated.sct)
AIM.integrated.sct <- RunUMAP(AIM.integrated.sct,reduction = "pca", dims = 1:30)
AIM.integrated.sct <- FindNeighbors(AIM.integrated.sct, reduction = "pca", dims = 1:30)
AIM.integrated.sct <- FindClusters(AIM.integrated.sct, resolution = 0.5)

p1 <- DimPlot(AIM.integrated.sct, group.by = "orig.ident")
p2 <-DimPlot(AIM.integrated.sct, label=TRUE, group.by  = "seurat_clusters")
p1   + p2 

p1 <- ggplot(AIM.integrated.sct@meta.data, aes(fill=factor(orig.ident), x=seurat_clusters,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") 
p1   + p2

saveRDS(AIM.integrated.sct, file = "202210223_AIM.integrated.sct.rds")
#AIM.integrated.sct <- readRDS(file ="~/S3data/Saved_Objects/202210223_AIM.integrated.sct.rds")

p1 <- RidgePlot(AIM.integrated.sct, features = c("CD3E","LYZ","MS4A1", "CD8A","GNLY","MIR155HG"), 
                ncol=4)
p1 <- RidgePlot(AIM.integrated.sct, features = c("IFNG","IFNA1", "IFNA2","MIR155HG"), 
                sort="decreasing")
p1   + p2
table(AIM.integrated.sct$seurat_clusters, AIM.integrated.sct$subject)
DefaultAssay(AIM.integrated.sct) <- "RNA"
FeaturePlot(object = AIM.integrated.sct, reduction = 'umap', features = c("MIR155HG"), cols = c("lightgray", "darkred")) 
FeaturePlot(object = AIM.integrated.sct, reduction = 'umap', features = c("IFNG", "TNF","IL2", "TNFRSF9", "TNFRSF4", "LTA", "CD200", "MIR155HG", "CD40LG", 
                                                                          "TNFSF8", "IL2RA", "TFRC"),
            cols = c("lightgray", "darkred"))                                                                                                              #AIM reactive genes 
FeaturePlot(object = AIM.integrated.sct, reduction = 'umap', features = c("CD3E", "MS4A1", 
                                                                          "GNLY", "LYZ", "NKG7", "CD4.1"), label = TRUE,
            cols = c("lightgray", "darkred"))   
DimPlot(AIM.integrated.sct, label = TRUE)
RidgePlot(AIM.integrated.sct, features = c("CD3E", "MS4A1", 
                                           "GNLY", "LYZ", "NKG7","CD8A", "CD4.1"), sort = "decreasing")
#-------------------------Read in integrated dataset--------------
AIM.integrated.sct <- readRDS(file ="~/S3data/Saved_Objects/202210223_AIM.integrated.sct.rds")
#DefaultAssay(AIM.integrated.sct) <- "ADT"
#AIM.integrated.sct <- NormalizeData(AIM.integrated.sct, normalization.method = 'CLR', margin = 2)

#---------------looking for batch effects in AIM integrated------------------------
#AIM.integrated.sct <- AIM.integrated.sct %>% PrepSCTFindMarkers() %>% as.SingleCellExperiment()
DefaultAssay(AIM.integrated.sct) <- "RNA"
AIM.integrated.sct$HTO_maxID <- as.factor(AIM.integrated.sct$HTO_maxID)
levels(AIM.integrated.sct$HTO_maxID)
AIM.integrated.sct$HTO_maxID <-  gsub('-7moPV', '-7moPostVax', AIM.integrated.sct$HTO_maxID)
AIM.integrated.sct$HTO_maxID <-  gsub('-1moPB', '-1moPostBoost', AIM.integrated.sct$HTO_maxID)
AIM.integrated.sct$HTO_maxID <-  gsub('Subject5', 'Subject5', AIM.integrated.sct$HTO_maxID)
AIM.integrated.sct$subject <-  gsub('Subject5', 'Subject5', AIM.integrated.sct$subject)
table(AIM.integrated.sct$subject,AIM.integrated.sct$timepoint)

AIM.integrated.sct@meta.data <- AIM.integrated.sct@meta.data %>%
  unite('subject_timepoint', subject,timepoint, remove=FALSE)
View(AIM.integrated.sct@meta.data)
AIM.integrated.sct.sce <- AIM.integrated.sct  %>% as.SingleCellExperiment()

## Determine the number of cells per sample
groups <- colData(AIM.integrated.sct.sce)[, c("subject_timepoint")]

AIM.integrated.sct.sce <- removeAltExps(AIM.integrated.sct.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk <- scuttle::aggregateAcrossCells(AIM.integrated.sct.sce, ids = colData(AIM.integrated.sct.sce)[, c("subject_timepoint")],use.assay.type = "counts")
class(pseudo_bulk)
dim(pseudo_bulk)

# transform, so rows are genes and columns are samples and make rownames as the sample IDs
pseudo_bulk <- assay(pseudo_bulk) 
#create metaData file 
metaData <- colnames(pseudo_bulk) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(endsWith(timepoint, "PostBreakThru") ~ "breakthrough" ,
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    grepl(  "Subject2", subject)  ~ "NOV",grepl(  "Subject3", subject)  ~ "JAN", grepl(  "Subject5", subject)  ~ "JAN", grepl(  "Subject6", subject)  ~ "JAN",
    grepl(  "Subject7", subject)  ~ "DEC",grepl(  "Subject8", subject)  ~ "MARCH",grepl(  "Subject2_PostBreakThru", subject)  ~ "MARCH",
    grepl(  "Subject9", subject) ~ "DEC", grepl(  "Subject12", subject) ~ "MARCH",grepl(  "Subject6_PostBreakThru", subject)  ~ "MARCH",
    grepl(  "Subject13", subject) ~ "JAN",grepl(  "Subject14", subject) ~ "MARCH", grepl(  "Subject1", subject)  ~ "SEPT",grepl(  "Subject4", subject)  ~ "SEPT", 
    grepl(  "Subject5_PostBreakThru", subject)  ~ "SEPT",grepl(  "Subject10", subject)  ~ "SEPT",grepl(  "Subject11", subject)  ~ "SEPT"))
metaData<- metaData %>%
  unite('timepoint_PriorCOVID', timepoint,PriorCOVID, remove=FALSE)
metaData

# Create DESeq2 object        
cluster_counts <- as.data.frame(as.matrix(pseudo_bulk[, which(colnames(pseudo_bulk) %in% metaData$subject_timepoint)])) 
#View(cluster_counts) #just making sure metaData and counts have same info
# Check that all of the row names of the metadata are the same 
#and in the same order as the column names of the counts in order to use as input to DESeq2
all(metaData$subject_timepoint == colnames(cluster_counts))         

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = metaData, 
                              design = ~ 1) 
#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)

z <- plotPCA(vsd, intgroup = "runDate")
z <- z + geom_point( size = 6)+ theme_classic()  
z
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")

z <- plotPCA(vsd, intgroup = "PriorCOVID")
z <- z + geom_point( size = 6)+ theme_classic()  
z
vsd$timepoint_PriorCOVID <- factor(vsd$timepoint_PriorCOVID, levels = c("7moPostVax_no","1moPostBoost_no","PostBreakThru_breakthrough","7moPostVax_yes","1moPostBoost_yes"))
z<- plotPCA(vsd, intgroup = "timepoint_PriorCOVID")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + 
  scale_color_manual(values = c(pre_vax, post_vax,   
                                breakthru, pre_infect, post_infect)) 
z
ggsave(z, filename = "./Images/Images/AllSample_PCA.pdf", width = 5.5, height = 5)
#no obvious outliers 

#-------------------------Subset on CD4 and CD8 T cells and cluster--------------
#subset on just CD4 and CD8 clusters 
AIM.integrated.T <-subset(AIM.integrated.sct, subset = seurat_clusters == "3" | 
                            seurat_clusters == "12"| seurat_clusters == "11"|
                            seurat_clusters == "13"|seurat_clusters == "2"|
                            seurat_clusters == "9"|seurat_clusters == "5"|
                            seurat_clusters == "15"|seurat_clusters == "22"|
                            seurat_clusters == "0"|
                            seurat_clusters == "16"|seurat_clusters == "7"|
                            seurat_clusters == "17")

DimPlot(AIM.integrated.T, label=TRUE, group.by  = "seurat_clusters")
VlnPlot(AIM.integrated.T,features = c("LYZ","MS4A1","CD19","CD3E","IL7R")) 

DefaultAssay(AIM.integrated.T)<-"integrated"
Assays(AIM.integrated.T)
#rm(AIM.integrated.sct) #free up memory in global environment

# These are now standard steps in the Seurat workflow for visualization and clustering
AIM.integrated.T <- RunPCA(AIM.integrated.T, verbose = TRUE)
ElbowPlot(AIM.integrated.T)
DimPlot(AIM.integrated.T, reduction = "pca")

AIM.integrated.T <- RunUMAP(AIM.integrated.T, dims = 1:5)
AIM.integrated.T <- FindNeighbors(AIM.integrated.T, reduction = "pca", dims = 1:10) #1:10 previously
AIM.integrated.T <- FindClusters(AIM.integrated.T, resolution = 0.5)  #0.5 previously

p1 <- DimPlot(AIM.integrated.T, split.by  = "sample", label = TRUE)
p2 <-DimPlot(AIM.integrated.T, group.by = "seurat_clusters", label=TRUE, label.size = 7) + NoLegend()
p1   + p2

p1 <- ggplot(AIM.integrated.T@meta.data, aes(fill=factor(PriorCOVID), x=seurat_clusters,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") 
p1   + p2
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.T_UMAP.pdf", width = 10, height = 8)
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.T_UMAP.jpeg", width = 10, height = 8)

table(AIM.integrated.T$seurat_clusters)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones

#saveRDS(AIM.integrated.T, file = "~/S3data/Saved_Objects/20221213_AIM.integrated.T.rds")
AIM.integrated.T <- readRDS(file = "~/S3data/20221213_AIM.integrated.T.rds")

DefaultAssay(AIM.integrated.T) <- "RNA"
AIM.integrated.T <- ScaleData(AIM.integrated.T)

pbmc.markers <- FindAllMarkers(AIM.integrated.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 2000 )
write.csv(pbmc.markers,"~/S3data/Export_csv/DEG_AIM.integrated.T_byCluster.csv", row.names = T)
#pbmc.markers <- read_csv("~/S3data/Export_csv/DEG_AIM.integrated.T_byCluster.csv" )

top3 <- pbmc.markers %>% subset(`p_val_adj` <= 0.01)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
View(top3)

DefaultAssay(AIM.integrated.T) <- 'RNA'
Idents(AIM.integrated.T) <- "seurat_clusters"
cluster.averages <- AverageExpression(AIM.integrated.T, return.seurat = TRUE)
z <-DoHeatmap(cluster.averages, features = top3$gene, draw.lines = FALSE) +NoLegend()
z
z+p2
ggsave(z, file = "~/S3data/Images/Images/AIMintegrated.T_Heatmap_of_clusters_CD4_CD8.pdf", width = 12, height = 12)

DimPlot(AIM.integrated.T, group.by = "seurat_clusters", label=TRUE, label.size = 7, split.by = "sample")
table(AIM.integrated.T$PriorCOVID,AIM.integrated.T$seurat_clusters)
DefaultAssay(AIM.integrated.T) <- "ADT"

Idents(AIM.integrated.T) <- "seurat_clusters"
integrated_cluster.IDs <- c("Naive CD4+ T", "CD4+ T","Naive CD8+ T", "CD8+ T",                                    #0-3
                            "Naive CD4+ T","CD4+ T", "AIM Reactive CD4+ T",  "CD8+ T",                         #4-7
                            "CD4+ T" ,"CD8+ T","Naive CD4+ T", "CD8+ T", "Treg",                        #8-12
                            "AIM Reactive CD8+ T", "Naive CD8+ T", "AIM Reactive CD4+ T")                           #13-15

names(integrated_cluster.IDs) <- levels(AIM.integrated.T)
AIM.integrated.T <- RenameIdents(AIM.integrated.T, integrated_cluster.IDs)
AIM.integrated.T@meta.data$integrated_cluster.IDs <- AIM.integrated.T@active.ident
DimPlot(AIM.integrated.T, label=F, group.by  = "integrated_cluster.IDs",raster=F)
DimPlot(AIM.integrated.T, label=F, group.by  = "integrated_cluster.IDs",raster=F, split.by = "sample")
table(AIM.integrated.T$integrated_cluster.IDs, AIM.integrated.T$subject)

cluster.averages <- AverageExpression(AIM.integrated.T, return.seurat = TRUE)
z <-DoHeatmap(cluster.averages, features = top3$gene, draw.lines = FALSE) +NoLegend()
z

#dump reads with MS4A1  
DefaultAssay(AIM.integrated.T) <- "RNA"
FeaturePlot(AIM.integrated.T,features = c("LYZ","MS4A1")) 
FeatureScatter(AIM.integrated.T, feature1 = "CD3E", feature2 = "MS4A1", group.by = "seurat_clusters")+ geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1)
#AIM.integrated.T.clean <- subset(AIM.integrated.T, subset = CD3E < 0.1 & MS4A1 > 0.1 |CD3E > 0.1 & MS4A1 < 0.1 | CD3E < 0.1 & MS4A1 < 0.1)
#rm(AIM.integrated.T.clean)
#RidgePlot(AIM.integrated.T, features = c("LYZ","MS4A1"), sort = "decreasing")

AIM.integrated.T.clean <- subset(AIM.integrated.T, subset = MS4A1 < 0.1 )#& LYZ < 0.1 )
rm(AIM.integrated.T)
#DimPlot(AIM.integrated.T.clean, label=TRUE, group.by  = "integrated_cluster.IDs", split.by = "sample")
#FeaturePlot(AIM.integrated.T.clean,features = c("LYZ","MS4A1")) 
#table(AIM.integrated.T.clean$integrated_cluster.IDs, AIM.integrated.T.clean$run)

#-------------------------Subset on CD4 T cells and redraw UMAP--------------
DefaultAssay(AIM.integrated.T.clean)<-"integrated"
#rm(AIM.integrated.sct) #free up memory in global environment

AIM.integrated.T.clean.CD4 <- subset(AIM.integrated.T.clean, integrated_cluster.IDs == "Naive CD4+ T" |
                                       integrated_cluster.IDs == "CD4+ T"|
                                       integrated_cluster.IDs == "Treg"|
                                       integrated_cluster.IDs == "AIM Reactive CD4+ T")

ElbowPlot(AIM.integrated.T.clean.CD4)
AIM.integrated.T.clean.CD4 <- RunUMAP(AIM.integrated.T.clean.CD4, dims = 1:8)

DimPlot(AIM.integrated.T.clean.CD4, group.by = "integrated_cluster.IDs")

Idents(AIM.integrated.T.clean.CD4) <- "integrated_cluster.IDs"
integrated_cluster.IDs <- c("Naive CD4+ T", "Nonnaive CD4+ T", "AIM Reactive CD4+ T",  "Treg")   

names(integrated_cluster.IDs) <- levels(AIM.integrated.T.clean.CD4)
AIM.integrated.T.clean.CD4 <- RenameIdents(AIM.integrated.T.clean.CD4, integrated_cluster.IDs)
AIM.integrated.T.clean.CD4@meta.data$integrated_cluster.IDs <- AIM.integrated.T.clean.CD4@active.ident
DimPlot(AIM.integrated.T.clean.CD4, label=T, group.by  = "integrated_cluster.IDs",raster=F)

pbmc.markers <- FindAllMarkers(AIM.integrated.T.clean.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 2000)
top3 <- pbmc.markers %>% subset(`p_val_adj` <= 0.01)
top3 <- top3 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top3$cluster <- factor(top3$cluster, levels = c("Naive CD4+ T", "Nonnaive CD4+ T", "Treg", "AIM Reactive CD4+ T"))
top3 <- arrange(top3, cluster)

DefaultAssay(AIM.integrated.T.clean.CD4) <- 'RNA'
Idents(AIM.integrated.T.clean.CD4) <- "integrated_cluster.IDs"
cluster.averages <- AverageExpression(AIM.integrated.T.clean.CD4, return.seurat = TRUE)
cluster.averages@active.ident <- factor(cluster.averages@active.ident, levels = c("Naive CD4+ T", "Nonnaive CD4+ T", "Treg", "AIM Reactive CD4+ T"))
z <-DoHeatmap(cluster.averages, features = top3$gene, draw.lines = FALSE, group.bar.height = 0.00,hjust = 0.08)+ scale_fill_viridis_c(direction = 1, option = "viridis") +
  theme(axis.text=element_text(colour="black",size=6)) 
z 
ggsave(z, filename = "./Images/Images/AIMintegrated.T_Heatmap_Top3Genes.pdf", width = 6, height = 8.5)

#-----------------------Adding contigs to metadata AIM.integrated.T--------------------------
setwd("~/S3data/Contig_annotations")
stim_NOV_TCR <- read.csv("20211114_ECCITEseq Pilot_Subject2/filtered_contig_annotations.csv")
stim_DEC_TCR <- read.csv("20211203_ECCITEseq_Subject9_Subject7/filtered_contig_annotations.csv")
stim1_JAN_TCR <- read.csv("20220131_ECCITEseq_AIM_sgg/stim1_outs_raw/filtered_contig_annotations.csv")
stim2_JAN_TCR <- read.csv("20220131_ECCITEseq_AIM_sgg/stim2_outs_raw/filtered_contig_annotations.csv")
stim1_MARCH_TCR <- read.csv("20220320_ECCITEseq_AIM/Stim1_outs_raw/filtered_contig_annotations.csv")
stim2_MARCH_TCR <- read.csv("20220320_ECCITEseq_AIM/Stim2_outs_raw/filtered_contig_annotations.csv")
stim3_MARCH_TCR <- read.csv("20220320_ECCITEseq_AIM/Stim3_outs_raw/filtered_contig_annotations.csv")
stim1_SEPT_TCR <- read.csv("20220927_ECCITEseq_AIM/Stim1_outs_raw/filtered_contig_annotations.csv")
stim2_SEPT_TCR <- read.csv("20220927_ECCITEseq_AIM/Stim2_outs_raw/filtered_contig_annotations.csv")
stim3_SEPT_TCR <- read.csv("20220927_ECCITEseq_AIM/Stim3_outs_raw/filtered_contig_annotations.csv")
stim4_SEPT_TCR <- read.csv("20220927_ECCITEseq_AIM/Stim4_outs_raw/filtered_contig_annotations.csv")
stim5_SEPT_TCR <- read.csv("20220927_ECCITEseq_AIM/Stim5_outs_raw/filtered_contig_annotations.csv")
# 
# stim_NOV_TCR$orig.ident <- "stim_NOV"
# stim_DEC_TCR$orig.ident <- "stim_DEC"
# stim1_JAN_TCR$orig.ident <- "stim1_JAN"
# stim2_JAN_TCR$orig.ident <- "stim2_JAN"
# stim1_MARCH_TCR$orig.ident <- "stim1_MARCH"
# stim2_MARCH_TCR$orig.ident <- "stim2_MARCH"
# stim3_MARCH_TCR$orig.ident <- "stim3_MARCH"
# stim1_SEPT_TCR$orig.ident <- "stim1_SEPT"
# stim2_SEPT_TCR$orig.ident <- "stim2_SEPT"
# stim3_SEPT_TCR$orig.ident <- "stim3_SEPT"
# stim4_SEPT_TCR$orig.ident <- "stim4_SEPT"
# stim5_SEPT_TCR$orig.ident <- "stim5_SEPT"

contig_list <- list(stim_NOV_TCR, stim_DEC_TCR, stim1_JAN_TCR, stim2_JAN_TCR,
                    stim1_MARCH_TCR,stim2_MARCH_TCR,stim3_MARCH_TCR,
                    stim1_SEPT_TCR,stim2_SEPT_TCR,stim3_SEPT_TCR,stim4_SEPT_TCR,stim5_SEPT_TCR)
View(contig_list)

combined <- combineTCR(contig_list, 
                       samples =  c("stim_NOV","stim_DEC","stim1_JAN","stim2_JAN","stim1_MARCH","stim2_MARCH","stim3_MARCH",
                                    "stim1_SEPT","stim2_SEPT","stim3_SEPT","stim4_SEPT","stim5_SEPT"),
                       ID = c("1", "2", "3", "4", "5", "6","7","8","9","10","11","12"),
                       cells ="T-AB")
View(contig_list)

for (i in seq_along(combined)) {
  combined[[i]] <- stripBarcode(combined[[i]], 
                                column = 1, connector = "_", num_connects = 4)
}

View(combined)

combined$stim_NOV_1$barcode_suffix <- "1"
combined$stim_DEC_2$barcode_suffix <- "2"
combined$stim1_JAN_3$barcode_suffix <- "3"
combined$stim2_JAN_4$barcode_suffix <- "4"
combined$stim1_MARCH_5$barcode_suffix <- "5"
combined$stim2_MARCH_6$barcode_suffix <- "6"
combined$stim3_MARCH_7$barcode_suffix <- "7"
combined$stim1_SEPT_8$barcode_suffix <- "8"
combined$stim2_SEPT_9$barcode_suffix <- "9"
combined$stim3_SEPT_10$barcode_suffix <- "10"
combined$stim4_SEPT_11$barcode_suffix <- "11"
combined$stim5_SEPT_12$barcode_suffix <- "12"

for (i in seq_along(combined)) {
  combined[[i]] <- combined[[i]]  %>% unite('barcode',barcode,barcode_suffix, remove=T)
}
View(combined)
saveRDS(combined, file = "~/S3data/Saved_Objects/combined.contigs.T.rds") 
#combined <- readRDS("~/S3data/Saved_Objects/combined.contigs.T.rds")

combined <- readRDS("~/S3data/Saved_Objects/combined.contigs.T.rds")

data_to_add <- rownames_to_column(AIM.integrated.T.clean.CD4@meta.data, "barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID")] %>% as.data.frame()
colnames(data_to_add) <- c("barcode","HTO_maxID")
for (i in seq_along(combined)) {
  combined[[i]] <- left_join(combined[[i]], data_to_add, by = "barcode")}

AIM.integrated.T.clean.CD4 <- combineExpression(combined, AIM.integrated.T.clean.CD4, cloneCall = "aa", chain = "TRB", group.by = "HTO_maxID", proportion = T,
                                                cloneTypes=c(Rare = 5e-4, Small = 0.001, Medium = 0.005, Large = 0.01, Hyperexpanded = 1))

slot(AIM.integrated.T.clean.CD4, "meta.data")$cloneType <- factor(slot(AIM.integrated.T.clean.CD4, "meta.data")$cloneType, 
                                                                  levels = c("Hyperexpanded (0.01 < X <= 1)",
                                                                             "Large (0.005 < X <= 0.01)",
                                                                             "Medium (0.001 < X <= 0.005)",
                                                                             "Small (5e-04 < X <= 0.001)",
                                                                             "Rare (0 < X <= 5e-04)",
                                                                             NA))
DimPlot(AIM.integrated.T.clean.CD4, group.by = "cloneType", raster = F, order = c("Hyperexpanded (0.01 < X <= 1)",
                                                                                  "Large (0.005 < X <= 0.01)",
                                                                                  "Medium (0.001 < X <= 0.005)",
                                                                                  "Small (5e-04 < X <= 0.001)",
                                                                                  "Rare (0 < X <= 5e-04)")) +
  scale_color_manual(values = colorblind_vector(5), na.value=NA) + 
  theme(plot.title = element_blank())

ggplot(AIM.integrated.T.clean.CD4@meta.data, aes(fill=factor(cloneType), x=integrated_cluster.IDs,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") + 
  scale_fill_manual(values = rev(colorblind_vector(5)), na.value="grey")


unique(AIM.integrated.T.clean.CD4$cloneType)

saveRDS(AIM.integrated.T.clean.CD4, file = "~/S3data/Saved_Objects/20221214_AIM.integrated.CD4T.rds") 
