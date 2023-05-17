library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
library(biovizBase)
library(data.table)
library(GenomicRanges)
library(tibble)
library(future)
#install.packages('qlcMatrix')
library(qlcMatrix)
#install.packages("reticulate")
library(reticulate)
library(EnsDb.Hsapiens.v86)
#BiocManager::install("Bioconductor/GenomeInfoDb", force = T)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install(version = "3.16")
library(BiocManager)
#BiocManager::install("scDblFinder")
library(scDblFinder)
library(glmGamPoi)

sessionInfo()
setwd("~/S3data/TEAseq")

#---------merge datasets------------ 
obj_stim1_1$dataset <- 'stim1_1'
obj_stim2_1$dataset <- 'stim2_1'
obj_stim3_1$dataset <- 'stim3_1'
obj_unstim_1$dataset <- 'unstim_1'

obj_stim1_2$dataset <- 'stim1_2'
obj_stim2_2$dataset <- 'stim2_2'
obj_stim3_2$dataset <- 'stim3_2'
obj_unstim_2$dataset <- 'unstim_2'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = obj_stim1_1,
  y = list(obj_stim2_1, obj_stim3_1, obj_unstim_1,
           obj_stim1_2,obj_stim2_2,obj_stim3_2,obj_unstim_2))

#View(combined@meta.data)
combined@meta.data <- combined@meta.data  %>%
  mutate(run = case_when(
    endsWith(dataset, "_1")  ~ "first",
    endsWith(dataset, "_2")  ~ "second" ))

combined@meta.data <- combined@meta.data  %>%
  mutate(condition = case_when(
    startsWith(dataset, "stim")  ~ "stim",
    startsWith(dataset, "unstim")  ~ "unstim"))

saveRDS(combined, file = "20230403_combined.no.Batch.Correction.RDS")

#-----------RNA integration-------------
#combined <- readRDS(file = "20230403_combined.no.Batch.Correction.RDS")

combined$dataset <- as.factor(combined$dataset)
levels(combined$dataset)
DefaultAssay(combined) <- "RNA"

# split the dataset into a list of the 4 seurat objects based on lane
combined <- SplitObject(combined, split.by = "dataset")

# normalize using SCTransform and identify variable features for each dataset independently
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
combined <- lapply(X = combined, FUN = SCTransform, method = "glmGamPoi")

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
#Run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = combined, nfeatures = 3000)
combined <- PrepSCTIntegration(object.list = combined, anchor.features = features)
combined <- lapply(X = combined, FUN = RunPCA, features = features)

#We then identify anchors using the FindIntegrationAnchors() function, which takes a list of 
#Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
immune.anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", 
                                         reference = c(1,4),  #reference datasets are form first run 
                                         anchor.features = features, dims = 1:30, 
                                         reduction = "rpca", k.anchor = 20)

combined.integrated.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
rm(immune.anchors,combined)
View(combined.integrated.sct)

combined.integrated.sct <- RunPCA(combined.integrated.sct)
ElbowPlot(combined.integrated.sct)

DefaultAssay(combined.integrated.sct) <- "integrated"
combined.integrated.sct <- RunUMAP(combined.integrated.sct, reduction = "pca", dims = 1:30, reduction.name = "gex_umap")
combined.integrated.sct <- FindNeighbors(combined.integrated.sct, reduction = "pca", dims = 1:30)
combined.integrated.sct <- FindClusters(combined.integrated.sct, resolution = 0.5, graph.name = "integrated_snn")
p1 <- DimPlot(combined.integrated.sct, group.by = "run", reduction = "gex_umap" )
p2 <-DimPlot(combined.integrated.sct, label=TRUE, group.by  = "seurat_clusters", reduction = "gex_umap")
p1   + p2 
DimPlot(combined.integrated.sct, label=TRUE, group.by  = "seurat_clusters", split.by = "run", reduction = "gex_umap")

p1 <- ggplot(combined.integrated.sct@meta.data, aes(fill=factor(run), x=seurat_clusters,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") 
p1   + p2
table(combined.integrated.sct$seurat_clusters)
saveRDS(combined.integrated.sct, file = "20230403_combined.rds") #ran up to here --- friday april 7 1:00AM 


rm(obj_stim1_1,obj_stim2_1, obj_stim3_1, obj_unstim_1,
   obj_stim1_2,obj_stim2_2,obj_stim3_2,obj_unstim_2)

#normalize ADT 
DefaultAssay(combined.integrated.sct) <- "ADT"
VariableFeatures(combined.integrated.sct) <- rownames(combined.integrated.sct[["ADT"]])
combined.integrated.sct <- NormalizeData(combined.integrated.sct, normalization.method = "CLR" , margin = 2)

DefaultAssay(combined.integrated.sct) <- "SCT"
FeaturePlot(combined.integrated.sct, features = c("adt_CD4", "CD8A","MS4A1","MIR155HG"), reduction = "gex_umap")
DimPlot(combined.integrated.sct, reduction = "gex_umap", label = T)
RidgePlot(combined.integrated.sct, features = c("adt_CD4")) #sort = "decreasing")
RidgePlot(combined.integrated.sct, features = c("MS4A1"), sort = "decreasing")

#subset on CD4 
CD4.integrated <- subset(combined.integrated.sct, seurat_clusters == 0| seurat_clusters == 1|
                          seurat_clusters == 2|seurat_clusters == 3|seurat_clusters == 4|
                          seurat_clusters == 6|seurat_clusters == 9|seurat_clusters == 11|
                          seurat_clusters == 12|seurat_clusters == 14|
                          seurat_clusters == 17|seurat_clusters == 19|seurat_clusters == 21|
                          seurat_clusters == 22)

#taking out clusters 5,7,8,10,13,16,15,18,20,23,24
DimPlot(CD4.integrated, reduction = "gex_umap",label = T)
FeaturePlot(CD4.integrated, features = c("adt_CD4"), reduction = "gex_umap", label = T)



#-----------ATAC integration-------------
# process the combined dataset
DefaultAssay(CD4.integrated) <- "ATAC"
CD4.integrated <- FindTopFeatures(CD4.integrated, min.cutoff = 10)
CD4.integrated <- RunTFIDF(CD4.integrated)
CD4.integrated <- RunSVD(CD4.integrated)
CD4.integrated <- RunUMAP(CD4.integrated, reduction = "lsi", dims = 2:10, reduction.name = "atac_umap")
View(CD4.integrated)
p1 <- DimPlot(CD4.integrated, group.by = "run", reduction = "atac_umap")
p1
ggsave(p1, filename = "~/S3data/Images/Images/TEAseq_CD4_atac_UMAP_not_integrated.jpeg", width = 5, height = 5)

#saveRDS(CD4.integrated, file = "20230407_CD4integrated.rds")
CD4.integrated <- readRDS(file = "20230407_CD4integrated.rds")

CD4.integrated.rlsi <- SplitObject(CD4.integrated, split.by = "dataset")
for (i in seq_along(CD4.integrated.rlsi)) {
  DefaultAssay(CD4.integrated.rlsi[[i]]) <- "ATAC"
  }

integration.anchors <- FindIntegrationAnchors(
  object.list = CD4.integrated.rlsi, 
  reference = c(1,4),
  anchor.features = rownames(CD4.integrated.rlsi$stim1_1),
  reduction = "rlsi",
  dims = 2:15
)

# integrate LSI embeddings
CD4.integrated.rlsi <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = CD4.integrated[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:15
)
rm(integration.anchors)
View(CD4.integrated.rlsi)
# create a new UMAP using the integrated embeddings
CD4.integrated.rlsi <- RunUMAP(CD4.integrated.rlsi, reduction = "integrated_lsi", dims = 2:15,
                                    reduction.name = "atac_umap" )
p2 <- DimPlot(CD4.integrated.rlsi, group.by = "condition", reduction = "atac_umap")
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))

DimPlot(CD4.integrated.rlsi, group.by = "condition", reduction = "atac_umap")

#saveRDS(CD4.integrated.rlsi, file = "20230407_CD4integrated.rlsi.rds")
CD4.integrated.rlsi <- readRDS(file = "20230407_CD4integrated.rlsi.rds")

#copy integrated LSI from duplicate seurat object to original object
CD4.integrated@reductions$integrated_lsi <- CD4.integrated.rlsi@reductions$integrated_lsi 

#overwrite redundant fragment paths 
DefaultAssay(CD4.integrated) <- "ATAC"
frags <- Fragments(CD4.integrated)
new.frags <- list(frags[[1]], frags[[10]],frags[[19]], frags[[28]], frags[[37]], frags[[46]],frags[[55]], frags[[64]])
View(new.frags)
Fragments(CD4.integrated) <- NULL
Fragments(CD4.integrated) <- new.frags
View(CD4.integrated)

CD4.integrated <- RunUMAP(CD4.integrated, reduction = "integrated_lsi", dims = 2:15,
                               reduction.name = "atac_integrated_umap" )
p2 <- DimPlot(CD4.integrated, group.by = "condition", reduction = "atac_integrated_umap")
p1 <- DimPlot(CD4.integrated, group.by = "condition", reduction = "CD4_umap")
z <- (p1 + ggtitle("RNA integrated")) | (p2 + ggtitle("ATAC integrated"))
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_atac_UMAP_integrates.jpeg", width = 10, height = 5)

p2 <- DimPlot(CD4.integrated, group.by = "seurat_clusters", reduction = "atac_integrated_umap")
p1 <- DimPlot(CD4.integrated, group.by = "seurat_clusters", reduction = "CD4_umap")
z <- (p1 + ggtitle("RNA integrated")) | (p2 + ggtitle("ATAC integrated"))
z
#----------gex clustering---------
DefaultAssay(CD4.integrated) <- "integrated"
CD4.integrated <- RunPCA(CD4.integrated)

ElbowPlot(CD4.integrated)

CD4.integrated <- RunUMAP(CD4.integrated, reduction = "pca", dims = 1:15, reduction.name = "gex_umap")
CD4.integrated <- FindNeighbors(CD4.integrated, reduction = "pca", dims = 1:15)
CD4.integrated <- FindClusters(CD4.integrated, resolution = 0.5, graph.name = "integrated_snn")
p1 <- DimPlot(CD4.integrated,  reduction = "gex_umap" )
p2 <-DimPlot(CD4.integrated, label=TRUE, reduction = "gex_umap")
z <- p1   + p2 
z
table(CD4.integrated$seurat_clusters)

DimPlot(CD4.integrated, label=TRUE, reduction = "gex_umap", group.by = "seurat_clusters")
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_UMAP.pdf", width = 10, height = 5)
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_UMAP.jpeg", width = 10, height = 5)

z<- DimPlot(CD4.integrated, label=TRUE, group.by  = "seurat_clusters", reduction = "gex_umap", split.by = "condition")
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_UMAP_Condition_Split.jpeg", width = 10, height = 5)

DefaultAssay(CD4.integrated) <- "SCT"
z<-FeaturePlot(CD4.integrated, features = c("MIR155HG", "TNFRSF9","IL2","IFNG","TNFRSF4","LTA"),order = T,
               reduction = "gex_umap", label = F)
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_UMAP_MIR155HG.jpeg", width = 5, height = 5)

ggplot(CD4.integrated@meta.data, aes(fill=factor(condition), x=seurat_clusters,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") 

table(CD4.integrated$condition)

VlnPlot(CD4.integrated, features = c("MIR155HG", "TNFRSF9", "FOXP3", "IFNG","GZMB"))
FeatureScatter(CD4.integrated, feature1 = "MIR155HG", feature2 = "TNFRSF9")

DimHeatmap(CD4.integrated, dims = 1:3, cells = 500, balanced = TRUE)
# Examine and visualize PCA results a few different ways
print(CD4.integrated[["pca"]], dims = 1:5, nfeatures = 5)

CD4.integrated.singlets <- subset(CD4.integrated, HTO_classification.global == "Singlet" |HTO_classification.global == "Negative" )
table(CD4.integrated.singlets$seurat_clusters, CD4.integrated.singlets$condition)

DefaultAssay(CD4.integrated.singlets) <- "RNA"
CD4.integrated.singlets <- CD4.integrated.singlets %>% NormalizeData() %>% ScaleData()
Idents(CD4.integrated.singlets) <- "seurat_clusters"
markers <- FindAllMarkers(CD4.integrated.singlets )
#write.csv(markers,"~/S3data/Export_csv/TEAseq_CD4_gex_umap_cluster_degs.csv", row.names = TRUE)
markers <- read.csv("~/S3data/Export_csv/TEAseq_CD4_gex_umap_cluster_degs.csv", row.names = 1)
markers.sig <- subset(markers, p_val_adj <= 0.05)

View(markers)
top5 <- markers.sig %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
CD4.integrated.avg <- AverageExpression(CD4.integrated.singlets, assay = "RNA", return.seurat = T)
DoHeatmap(CD4.integrated.avg, features = top5$gene, draw.lines = F, size = 5) +   
  theme(text =  element_text(size = 8))+ 
  scale_fill_viridis_c(direction = 1, option = "viridis")
scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = F)
    
table(CD4.integrated.singlets$restricted_clustering, CD4.integrated.singlets$HTO_maxID, CD4.integrated.singlets$condition )

Idents(CD4.integrated.singlets) <- "seurat_clusters"
# add annotations
CD4.integrated.singlets <- RenameIdents(CD4.integrated.singlets, '0' = 'Naive CD4+ T',
                                        '2' = 'Naive CD4+ T', '3' = 'Naive CD4+ T',
                                        '6' = 'Naive CD4+ T',
                                        '12' = 'Naive CD4+ T','18' = 'Naive CD4+ T')
CD4.integrated.singlets <- RenameIdents(CD4.integrated.singlets, '1' = 'Nonnaive CD4+ T','4' = 'Nonnaive CD4+ T','5' = 'Nonnaive CD4+ T',
                                        '7' = 'Nonnaive CD4+ T',
                                        '8' = 'Nonnaive CD4+ T','9' = 'Nonnaive CD4+ T','10' = 'Nonnaive CD4+ T','11' = 'Nonnaive CD4+ T',
                                        '14' = 'Nonnaive CD4+ T','16' = 'Nonnaive CD4+ T')
CD4.integrated.singlets <- RenameIdents(CD4.integrated.singlets, '13' = 'Treg', '15' = 'Treg')
CD4.integrated.singlets <- RenameIdents(CD4.integrated.singlets, '17' ='AIM Reactive CD4+ T')

table(Idents(CD4.integrated.singlets), CD4.integrated.singlets$condition)

Idents(CD4.integrated.singlets) <- factor(Idents(CD4.integrated.singlets), levels = c("Naive CD4+ T","Nonnaive CD4+ T",
                                                                                      "AIM Reactive CD4+ T","Treg"))

z <- DimPlot(CD4.integrated.singlets, label = F, repel = TRUE, reduction = "gex_umap")
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_UMAP_named.jpeg", width = 7, height = 5)

CD4.integrated.singlets.subset <- subset(x = CD4.integrated.singlets, downsample = 1000)
table(Idents(CD4.integrated.singlets.subset))

DefaultAssay(CD4.integrated.singlets.subset) <- "ATAC"
idents.plot <- Idents(CD4.integrated.singlets.subset)
p1 <- CoveragePlot(
  object = CD4.integrated.singlets.subset,
  region = "MIR155HG",extend.upstream = 10000, extend.downstream = 10000,
  expression.assay = "RNA",idents = idents.plot) 

p2 <- CoveragePlot(
  object = CD4.integrated.singlets.subset,
  region = "LTA", 
  expression.assay = "RNA", idents = idents.plot)

z <-patchwork::wrap_plots(p1, p2 ,ncol = 1) & scale_color_manual(values = c("darkgrey","grey")) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_coverage_plots.jpeg", width = 6, height = 6)

DefaultAssay(CD4.integrated.singlets.subset)<- "RNA"
VlnPlot(CD4.integrated.singlets.subset, features=c("TNFRSF9","MIR155HG"), pt.size = 1)


markers <- FindAllMarkers(CD4.integrated.singlets)
write.csv(markers,"~/S3data/Export_csv/TEAseq_CD4_gex_umap_integrated_clusterIDs_degs.csv", row.names = TRUE)
markers.sig <- subset(markers, p_val_adj <= 0.001)
View(markers.sig)
top5 <- markers.sig %>% group_by(cluster) %>% top_n(n = 35, wt = avg_log2FC)
CD4.integrated.avg <- AverageExpression(CD4.integrated.singlets, assay = "RNA", return.seurat = T)
z <- DoHeatmap(CD4.integrated.avg, features = top5$gene, draw.lines = F, size = 5) +   
  theme(text =  element_text(size = 4))+ 
  scale_fill_viridis_c(direction = 1, option = "viridis")
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_HeatMap_topgenes.pdf", width = 5, height = 7)

DefaultAssay(CD4.integrated.singlets) <- "RNA"
z<- FeaturePlot(CD4.integrated.singlets, reduction = "gex_umap", features= c("MIR155HG","TNFRSF9","TNFRSF4", "LTA"), ncol = 2, cols = 
                  c("lightgrey","darkred"), order = T, raster = F) & NoLegend() & NoAxes()
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_CD4_gex_UMAP_featurePlot.jpeg", width = 6, height = 5)
DimPlot(CD4.integrated.singlets, label=T)
table(CD4.integrated.singlets$seurat_clusters, CD4.integrated.singlets$condition)
#-------Reading in leading edge, hallmark gene sets, and vaccine and infection signatures------
#leading edge of G2M and mitotic spindle / IFNa and IFNg 
MitoticS_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_MITOTIC_SPINDLE.tsv", sep="\t")
MitoticS_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_MITOTIC_SPINDLE.tsv", sep="\t")
G2M_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_G2M_CHECKPOINT.tsv", sep="\t")
G2M_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_G2M_CHECKPOINT.tsv", sep="\t")
IFNa_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")
IFNa_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")
IFNg_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")
IFNg_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")

prolif_MS.pre <- subset(MitoticS_pre,CORE.ENRICHMENT == "Yes")
prolif_MS.post <- subset(MitoticS_post,CORE.ENRICHMENT == "Yes")
prolif_G2M.pre <- subset(G2M_pre,CORE.ENRICHMENT == "Yes")
prolif_G2M.post <- subset(G2M_post,CORE.ENRICHMENT == "Yes")
prolif_leading_edge <- c(prolif_MS.pre$SYMBOL, prolif_MS.post$SYMBOL,prolif_G2M.pre$SYMBOL,prolif_G2M.post$SYMBOL)
prolif_leading_edge <- unique(prolif_leading_edge)

ISGa.pre <- subset(IFNa_pre,CORE.ENRICHMENT == "Yes")
ISGa.post <- subset(IFNa_post,CORE.ENRICHMENT == "Yes")
ISGg.pre <- subset(IFNg_pre,CORE.ENRICHMENT == "Yes")
ISGg.post <- subset(IFNg_post,CORE.ENRICHMENT == "Yes")
ISG_leading_edge <- c(ISGa.pre$SYMBOL, ISGa.post$SYMBOL,ISGg.pre$SYMBOL, ISGg.post$SYMBOL)
ISG_leading_edge <- unique(ISG_leading_edge)

#full hallmark gene set
hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "H" # Only hallmark gene sets
)

hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)

ISGa <- hallmarks_list$HALLMARK_INTERFERON_ALPHA_RESPONSE
ISGg <-  hallmarks_list$HALLMARK_INTERFERON_GAMMA_RESPONSE
ISG <- c(ISGa, ISGg)
ISG <- unique(ISG)

prolif_MS<- hallmarks_list$HALLMARK_MITOTIC_SPINDLE
prolif_G2M <- hallmarks_list$HALLMARK_G2M_CHECKPOINT
prolif <- c(prolif_MS, prolif_G2M)
prolif <- unique(prolif)


INFLAM <- hallmarks_list$HALLMARK_INFLAMMATORY_RESPONSE
infect_hallmark_sig <- c(ISGa, ISGg,INFLAM)
infect_hallmark_sig <- unique(infect_hallmark_sig)
write.table(infect_hallmark_sig, file = "~/S3data/Gene_sets/infect_hallmark_sig.txt",row.names = FALSE, 
            col.names = F, quote =F,sep = "\t")

NFkb <- hallmarks_list$HALLMARK_TNFA_SIGNALING_VIA_NFKB
vax_hallmark_sig <- c(prolif_MS, prolif_G2M,NFkb)
vax_hallmark_sig <- unique(vax_hallmark_sig)
write.table(vax_hallmark_sig, file = "~/S3data/Gene_sets/vax_hallmark_sig.txt",row.names = FALSE, 
            col.names = F, quote =F,sep = "\t")

PriorCOVID.CD4 <- read.csv(file = "~/S3data/Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost.DESeq.csv", row.names = 1)
Infection_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Vaccine_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )

#-------subset on Spike-specific cluster----------------
CD4.integrated.singlets@meta.data <- CD4.integrated.singlets@meta.data  %>%
  mutate(timepoint = case_when(
    endsWith(HTO_maxID, "tion")  ~ "PostInfection",
    endsWith(HTO_maxID, "Dose")  ~ "PostVaccination"
  ))
Spike_Specific_CD4 <- subset ( CD4.integrated.singlets, seurat_clusters == 17 & condition == "stim")
table(Idents(Spike_Specific_CD4))
table(Spike_Specific_CD4$timepoint)

#----------DEG BETWEEN nonnaive CD4---------
CD4.integrated.singlets$integrated_clusterIDs <-  CD4.integrated.singlets@active.ident
nonaiveCD4 <- subset (CD4.integrated.singlets, integrated_clusterIDs == "Nonnaive CD4+ T" & condition == "unstim")
#View(nonaiveCD4)
Idents(nonaiveCD4)<- "timepoint"
DefaultAssay(nonaiveCD4) <- "RNA"

PriorCOVID.CD4 <- FindMarkers(nonaiveCD4,   ident.1 = "PostInfection",
                              ident.2 = "PostVaccination", test.use="wilcox") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.15 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.15 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"
#View(PriorCOVID.CD4)
Experienced.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Naive.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
Experienced_UP_Naive_DOWN_CD4.DESeq <- PriorCOVID.CD4
top <- rbind(rownames(Experienced.AIMreactive.CD4.DEseq),rownames(Naive.AIMreactive.CD4.DEseq))
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ #scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(vax,"lightgrey",infect))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+
  NoLegend()+labs(y = "-log10 p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z 


#---------find differential gene expression--------------
Idents(Spike_Specific_CD4)<- "timepoint"
DefaultAssay(Spike_Specific_CD4) <- "RNA"

PriorCOVID.CD4 <- FindMarkers(Spike_Specific_CD4,   ident.1 = "PostInfection",
                              ident.2 = "PostVaccination", test.use="DESeq2") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.15 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.15 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"
#View(PriorCOVID.CD4)
x <- rownames_to_column(PriorCOVID.CD4, "genes")
View(x)
Experienced.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Naive.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
Experienced_UP_Naive_DOWN_CD4.DESeq <- PriorCOVID.CD4
#View(Experienced.AIMreactive.CD4.DEseq)
top <- c("MTRNR2L8","IFI44L","HLA-C","DDX60","DDX60L")
#top <- c("HLA-B","CCL5","GZMH","NKG7","IFITM1","CCL4","IFITM3","CCL3","GZMB","IFI6","LTB","CST7",
#         "CCR7","DDIT4","NFKB2")
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ #scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(vax,"lightgrey",infect))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+
  NoLegend()+labs(y = "-log10 p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z 
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_SpikeSpecific_Volcano.pdf", width = 8, height =7)



infection_list <- list(ISG,  ISG_leading_edge, prolif, prolif_leading_edge,
                       rownames(Infection_imprint), rownames(Vaccine_imprint),
                       rownames(Experienced.AIMreactive.CD4.DEseq), rownames(Naive.AIMreactive.CD4.DEseq))
names(infection_list) = c("ISG","ISG_leading_edge",
                          "Prolif","prolif_leading_edge",
                          "Infection_imprint", "Vaccine_imprint", 
                          "infection_gex_up","vaccine_gex_up")

ggvenn(infection_list, c("ISG_leading_edge","infection_gex_up"), auto_scale = T)            # Pairwise venn diagram
ggvenn(infection_list, c("Prolif","vaccine_gex_up"), auto_scale = T) 

Experienced.AIMreactive.CD4.DEseq <-  Experienced.AIMreactive.CD4.DEseq %>% 
  add_column(overlap_ISG =
               rownames(Experienced.AIMreactive.CD4.DEseq) %in% ISG)

DefaultAssay(Spike_Specific_CD4) <- "RNA"
z <- VlnPlot(Spike_Specific_CD4, features = c("IFI44L"),ncol = 1,
             group.by  = "timepoint", pt.size = 0.0, y.max = 5) & 
  scale_fill_manual(values = c(infect, vax))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "~/S3data/Images/Images/TEA_seq_VlnPlot_DEGs_IFI44L.pdf", width = 2.5, height = 2)

DefaultAssay(Spike_Specific_CD4) <- "RNA"
z <- VlnPlot(Spike_Specific_CD4, features = c("HLA-C"),ncol = 1,
             group.by  = "timepoint", pt.size = 0.0, y.max = 5) & 
  scale_fill_manual(values = c(infect, vax))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "~/S3data/Images/Images/TEA_seq_VlnPlot_DEGs_HLA-C.pdf", width = 2.5, height = 2)

#-------pseduo bulk of Spike-specific GEX---------
DefaultAssay(Spike_Specific_CD4) <- "RNA"
AIM.reactive.CD4.sce <- Spike_Specific_CD4 %>% as.SingleCellExperiment() #using RNA counts !!!
table(Spike_Specific_CD4$HTO_maxID, Spike_Specific_CD4$run)
## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$HTO_maxID)
groups <- colData(AIM.reactive.CD4.sce)[, c("HTO_maxID")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("HTO_maxID")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "-",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject5")  ~ "first",endsWith(subject, "Subject6")  ~ "first", endsWith(subject, "Subject8")  ~ "first", endsWith(subject, "Subject9")  ~ "first",
    endsWith(subject, "Subject1")  ~ "second",endsWith(subject, "Subject24")  ~ "second",endsWith(subject, "Subject17")  ~ "second",
    endsWith(subject, "Subject7") ~ "second", endsWith(subject, "Subject23") ~ "second",  endsWith(subject, "Subject20") ~ "second", endsWith(subject, "Subject12") ~ "second"))
metaData

# Create DESeq2 object        
pseudo_bulk_CD4 <- assay(pseudo_bulk_CD4) 
cluster_counts <- as.data.frame(as.matrix(pseudo_bulk_CD4[, which(colnames(pseudo_bulk_CD4) %in% metaData$subject_timepoint)])) 
#View(cluster_counts) #just making sure metaData and counts have same info
# Check that all of the row names of the metadata are the same 
#and in the same order as the column names of the counts in order to use as input to DESeq2
all(metaData$subject_timepoint == colnames(cluster_counts))         

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = metaData, 
                              design = ~ timepoint) 

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "timepoint")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(infect,vax)) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~timepoint, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "subject")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "timepoint")+ geom_point( size = 2)+ theme_classic() +  scale_color_manual(values = c(infect,vax)) 
z<- plotPCA(vsd, intgroup = "timepoint")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + 
  scale_color_manual(values = c(infect,vax))  #"#ffa527","#56e39d", "#716be4" )) 
# scale_x_continuous(limits = c(-25,25)) + scale_y_continuous(limits = c(-10,10))
z


#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 4 # filter for genes with at least 10 counts in 3/12 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("timepoint","PostInfection","PostSecondDose"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/TEAseq_PostInfection_UP_PostSecondDose_DOWN_CD4_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")]
write.table(geneRank, file = "~/S3data/Gene_sets/TEAseq_PostInfection_UP_PostSecondDose_DOWN_CD4_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F,sep = "\t")

#--------GSEA results------------
PostInfect_UP_PostVax_DOWN_pos <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/TEAseq_PostInfection_UP_PostSecondDose_DOWN_PseudoBulk_DESeq/gsea_report_for_na_pos_1682451607584.tsv")
PostInfect_UP_PostVax_DOWN_neg <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/TEAseq_PostInfection_UP_PostSecondDose_DOWN_PseudoBulk_DESeq/gsea_report_for_na_neg_1682451607584.tsv")
GSEA <- rbind(PostInfect_UP_PostVax_DOWN_pos,PostInfect_UP_PostVax_DOWN_neg)
GSEA$NAME<-gsub("HALLMARK_","",as.character(GSEA$NAME))
GSEA$NAME<- str_to_title(GSEA$NAME) 
View(GSEA)

plotGSEAinfection <- function (Post,postPathway, pathwayName, title, leftLabel, rightLabel, cohortCompare, legendUpperRight)
{
  if (cohortCompare == F)     {     colorGrob1 <- "black"; colorGrob2 <- "black"     }
  if (cohortCompare == T)   { colorGrob1 <- "#ffa527"; colorGrob2 <-  "black"}
  annotationInfo <- paste0("NES: ", round(Post$NES[grep(pathwayName,Post$NAME)],2), "\n", "FDR: ", formatC(Post$`FDR q-val`[grep(pathwayName,Post$NAME)], format = "e", digits = 1))
  if (legendUpperRight == T)   {   grob1 = grobTree(textGrob(annotationInfo,  x=0.68,  y=0.82, hjust=0, gp=gpar(col= "black", fontsize=20))) }
  if (legendUpperRight == F)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.25, hjust=0, gp=gpar(col="black", fontsize=24)))     }
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.97,hjust=0, gp=gpar(col=colorGrob1, fontsize=20)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.7,y=0.97,hjust=0, gp=gpar(col=colorGrob2, fontsize=20))) 
  return(
    ggplot(data=postPathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES) ) + geom_line(color= "black", size=1) + 
      geom_rug(sides="b", size=0.75, alpha=0.5, color= "black" ) + theme_bw() +
      ggtitle(title) + ylab("Enrichment score") + xlab("Rank in gene list") + 
      theme(axis.text = element_text(size=20,hjust = 0.75, color="black"), axis.title = element_text(size=24,hjust = 0.5, color="black"), plot.title = element_text(size=32,hjust = 0.5))+
      annotation_custom(grob1) + geom_hline(yintercept = 0) + 
      annotation_custom(left_grob) + annotation_custom(right_grob) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
      geom_hline(yintercept = 0, size=0.5)
  )
}

IFNa <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/TEAseq_PostInfection_UP_PostSecondDose_DOWN_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")

a<- plotGSEAinfection(GSEA, IFNa, pathwayName = "Interferon_alpha_response", title = "Interferon Alpha Response", 
                      leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
a

IFNg <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/TEAseq_PostInfection_UP_PostSecondDose_DOWN_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")
b<- plotGSEAinfection(GSEA,IFNg, pathwayName = "Interferon_gamma_response", title = "Interferon Gamma Response", 
                      leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
b

ggsave(a, filename = "./Images/Images/GSEA_IFNa_PostInfect.pdf", width = 7, height = 4.5)
ggsave(b, filename = "./Images/Images/GSEA_IFNg_PostInfect.pdf", width = 7, height = 4.5)


#-----------GSVA infection and vaccine imprints------------
DefaultAssay(Spike_Specific_CD4) <- "RNA"
#AIM.reactive.CD4.preBoost <- NormalizeData(AIM.reactive.CD4.preBoost)
ActivatedT.mxt.RNA <- as.data.frame(Spike_Specific_CD4@assays$RNA@scale.data) #36601 genes
#View(ActivatedT.mxt)
ActivatedT.mxt.RNA <- as.matrix(ActivatedT.mxt.RNA)
#View(ActivatedT.mxt.RNA)

imprints <- list(rownames(Infection_imprint),rownames(Vaccine_imprint))
names(imprints) = c("Infection_imprint","Vaccine_imprint")
View(imprints)

prolif_ISG <- list(prolif, ISG)
names(prolif_ISG) = c("prolif","ISG")
View(prolif_ISG)

hallmark_sig <- list(infect_hallmark_sig, vax_hallmark_sig)
names(hallmark_sig) = c("infect_hallmark_sig","vax_hallmark_sig")
View(hallmark_sig)

gsva.es.all <- gsva(ActivatedT.mxt.RNA, hallmark_sig, verbose=TRUE, method="gsva")

dim(gsva.es.all)

metaData <- colnames(ActivatedT.mxt.RNA)
data_to_add <- rownames_to_column(Spike_Specific_CD4@meta.data ,"barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID","timepoint")]
colnames(data_to_add) <- c("barcodes","HTO_maxID","timepoint") 

metaData <- merge(metaData, data_to_add, by.x = 1, by.y = 1)
## Determine the number of cells per sample

colnames(metaData) <- c("barcode", "HTO_maxID","timepoint") 
metaData <- metaData %>% arrange(HTO_maxID)

gsva.es.all <- gsva.es.all[ , metaData$barcode]

metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("barcode")

gsva.results <- t(gsva.es.all)
data <- merge(metaData, gsva.results, by.x = 0, by.y = 0) %>% as.data.frame()

z <- ggplot(data, aes(x=infect_hallmark_sig,y=vax_hallmark_sig, color = timepoint)) +
  theme_classic() + geom_point(size=1) + geom_density2d(binwidth = 0.25,color="grey", size = 0.25) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(infect,vax)) + labs(x =  "Infection signature enrichment score",
                                                                                           y="Vaccine signature enrichment score") +  
  facet_grid( ~ timepoint) + NoLegend() + theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z
#--------Call peaks--------
#saveRDS(Spike_Specific_CD4, file = "Spike_Specific_CD4_TEAseq.rds" )  #SAVED WITH NEW PEAKS IN "peaks" ASSAY 04/12/23
#Spike_Specific_CD4 <- readRDS(file = "~/S3data/TEAseq/Spike_Specific_CD4_TEAseq.rds" )
Spike_Specific_CD4 <- readRDS("~/S3data/TEAseq/20230419_Spike_Specific_CD4_TEAseq.rds")


#rewrite Spike_Specific_CD4@assays$ATAC
#call peaks on subsetted cluster 
DefaultAssay(Spike_Specific_CD4) <- "ATAC"
table(Spike_Specific_CD4$orig.ident)
peaks <- CallPeaks(object = Spike_Specific_CD4, macs2.path = "/usr/local/bin/macs2") # started at 8:17 PM 
       
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
           
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(Spike_Specific_CD4),  
  features = peaks,   
  cells = colnames(Spike_Specific_CD4)
)         
 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"   

# create a new assay using the MACS2 peak set and add it to the Seurat object
Spike_Specific_CD4[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(Spike_Specific_CD4),
  annotation = annotation
)
Annotation(Spike_Specific_CD4[["peaks"]])

#saveRDS(Spike_Specific_CD4, file = "20230419_Spike_Specific_CD4_TEAseq.rds" ) 

DefaultAssay(Spike_Specific_CD4) <- "peaks"
Spike_Specific_CD4 <- FindTopFeatures(Spike_Specific_CD4, min.cutoff = 5)
Spike_Specific_CD4 <- RunTFIDF(Spike_Specific_CD4)
Spike_Specific_CD4 <- RunSVD(Spike_Specific_CD4)

Spike_Specific_CD4 <- RunUMAP(Spike_Specific_CD4, reduction = "lsi", dims = 3:10, reduction.name = "atac_umap")
View(Spike_Specific_CD4)
p1 <- DimPlot(Spike_Specific_CD4, group.by = "timepoint", reduction = "atac_umap")
p1

VlnPlot(Spike_Specific_CD4, features = c("LSI_1","LSI_2","LSI_3","LSI_4","LSI_5",
                                         "LSI_6","LSI_7","LSI_8","LSI_9","LSI_10"), group.by = "run")

#run differential accesibility analysis
Idents(Spike_Specific_CD4) <- (Spike_Specific_CD4$timepoint)
table(Spike_Specific_CD4$timepoint)
DefaultAssay(Spike_Specific_CD4) <- 'peaks'

da_peaks <- FindMarkers(
  object = Spike_Specific_CD4,
  ident.1 = "PostInfection",
  ident.2 = "PostVaccination",
  test.use = 'LR', min.pct = 0.05, 
  latent.vars = "nCount_peaks"
) 

#View(da_peaks)
da_peaks_sig <- subset(da_peaks, da_peaks$p_val < 0.05)
da_peaks_sig$diffexpressed <- "NO"
da_peaks_sig$diffexpressed[da_peaks_sig$avg_log2FC > 0 & 
                             da_peaks_sig$p_val < 0.05] <- "UP"
da_peaks_sig$diffexpressed[da_peaks_sig$avg_log2FC < 0 & 
                             da_peaks_sig$p_val < 0.05] <- "DOWN"

da_peaks_rownames <- rownames(da_peaks_sig)
closest_genes <- ClosestFeature(Spike_Specific_CD4, regions = da_peaks_rownames)
#View(da_peaks_sig)

da_peaks_sig <- rownames_to_column(da_peaks_sig, "query_region")
annotated_peaks <- merge(da_peaks_sig, closest_genes, by ="query_region")
#View(annotated_peaks)

#
infect_open <- subset(annotated_peaks,avg_log2FC > 0 )
vax_open <- subset(annotated_peaks,avg_log2FC < 0 )

#install.packages("ggvenn")
library(ggvenn)

infection_list <- list(ISG,  ISG_leading_edge, infect_open$gene_name, prolif, prolif_leading_edge, vax_open$gene_name, rownames(Infection_imprint), rownames(Vaccine_imprint))
names(infection_list) = c("ISG","ISG_leading_edge","infect_open", 
                          "Prolif","prolif_leading_edge","Vax_open", "Infection_imprint", "Vaccine_imprint")
#View(infection_list)
ggvenn(infection_list, c("Infection_imprint","infect_open"), auto_scale	=T)            # Pairwise venn diagram

infect_open <-  infect_open%>% 
  add_column(overlap_ISG =
               infect_open$gene_name %in% ISG)
infect_open <-  infect_open%>% 
  add_column(overlap_Infection_sign =
               infect_open$gene_name %in% rownames(Infection_imprint))
View(infect_open)


vax_open <-  vax_open%>% 
  add_column(overlap_prolif =
               vax_open$gene_name %in% prolif)
vax_open <-  vax_open%>% 
  add_column(overlap_vax_sig =
               vax_open$gene_name %in% rownames(Vaccine_imprint))
vax_open <-  vax_open%>% 
  add_column(overlap_vax_hallmark_sig =
               vax_open$gene_name %in% vax_hallmark_sig)
View(vax_open)

vax <- "#ffc471"
infect <- "#b5b2f1"

z<- ggplot(da_peaks_sig, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(infect,vax))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "# of differentially accessible peaks", x ="") + 
  geom_text(stat='count',aes(label=..count..), vjust = -0.5, size = 10) +  scale_y_continuous(limits = c(0,1500))+NoLegend()
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_summary_BarGraph.jpeg.jpeg", width = 4, height =6)

annotated_peaks$type <- factor(annotated_peaks$type, levels = c("gap","cds","exon","utr"))
z<- ggplot(annotated_peaks, aes(x=type))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(infect,vax))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20),
                        axis.title.x =element_text(colour="black",size=20)) +labs(y = "# of differentially accessible peaks", x ="Category") + 
  scale_y_continuous(limits = c(0,1000))+NoLegend()
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_categories_BarGraph.jpeg", width = 6, height =6)

#-------Visualizing DARs----------
DefaultAssay(Spike_Specific_CD4) <- "peaks"
# first compute the GC content for each peak
Spike_Specific_CD4 <- RegionStats(Spike_Specific_CD4, genome = BSgenome.Hsapiens.UCSC.hg38)
View(vax_open)


#STAT4
ranges.show <- StringToGRanges(c("chr2-191109915-191110190"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "chr2-191108915-191111190", 
    links = F, 
    extend.upstream = 10000,
    extend.downstream = 5000 ) 
  z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
  z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_STAT4.pdf", width = 8, height =2.5)

# CDK6
ranges.show <- StringToGRanges(c("chr7-92648127-92648429","chr7-92731162-92732911"))
ranges.show$color <- "lightblue"
z<- CoveragePlot(
      object = Spike_Specific_CD4, region.highlight	= ranges.show,
      region = "chr7-92647127-92733911",  features = "CDK6",  expression.assay = "RNA",
      links = F, 
      extend.upstream = 100,
      extend.downstream = 1000) 
z<- z & scale_fill_manual(values = c(vax,infect))   & scale_color_manual(values = c("darkgrey","grey")) 
z
 
#CCR7 
ranges.show <- StringToGRanges(c("chr17-4Subject123465-4Subject125867"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                region = "CCR7",  features = "CCR7",expression.assay = "RNA",
                 links = F, 
                   extend.upstream = 0,
                   extend.downstream = 5000) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#GRAMD1A
ranges.show <- StringToGRanges(c("chr19-34999474-35001405"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                   region = "GRAMD1A",  features = "GRAMD1A",
                   links = F, annotation = T,
                   extend.upstream = 0,
                   extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#DDIT4
ranges.show <- StringToGRanges(c("chr10-72259985-72261851"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "DDIT4",  features = "DDIT4",
    links = F, 
    extend.upstream = 20000,
    extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#infect_open 
View(infect_open)
#LY6E
ranges.show <- StringToGRanges(c("chr8-143023888-143024658"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "chr8-143023888-143024658",
    links = F, 
    extend.upstream = 7000,
    extend.downstream = 1500) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_LY6E.pdf", width = 8, height =2.5)

#FPR1
ranges.show <- StringToGRanges(c("chr19-51784286-51785431"))
ranges.show$color <- "darkgrey"
  DefaultAssay(Spike_Specific_CD4) <- "peaks"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, #region.highlight	= ranges.show,
    region = "chr19-51784286-51785431",  
    links = F, annotation = T,peaks = F,
    extend.upstream = 3000,
    extend.downstream = 6000) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey"))
z 
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_FPR1.pdf", width = 8, height =2.5)


#IFIT1
ranges.show <- StringToGRanges(c("chr10-89391631-89393326"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, #region.highlight	= ranges.show,
    region = "IFIT1", 
    links = F, 
    extend.upstream = 1000,
    extend.downstream = 1000) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_IFIT1.pdf", width = 8, height =2.5)

#LGALS3BP
ranges.show <- StringToGRanges("chr17-78970987-78971608")
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                   region = "LGALS3BP",
                   links = F, 
                   extend.upstream = 0,
                   extend.downstream = 0) 
  z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
  z
  ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_LGALS3BP.pdf", width = 8, height =2.5)

#SAMD9L
ranges.show <- StringToGRanges("chr7-93165096-93165863")
ranges.show$color <- "darkgrey"
z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                     region = "chr7-93165096-93165863",
                     links = F, 
                     extend.upstream = 19000,
                     extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_DAR_SAMD9L.pdf", width = 8, height =2.5)
    
#IFIT1 and IFIT3
DefaultAssay(Spike_Specific_CD4) <- "peaks"
ranges.show <- StringToGRanges(c("chr10-89392023-89392963", "chr10-89327186-89328396"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "chr10-89327186-89392963",
    links = F, peaks = F,
    extend.upstream = 100000,
    extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z


#B2M
ranges.show <- StringToGRanges(c("chr15-44712429-44713192"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "B2M",  features = "B2M",
    links = F, peaks = F,
    extend.upstream = 0,
    extend.downstream = 0) 
  z<- z & scale_fill_manual(values = c(infect,vax)) & scale_color_manual(values = c("darkgrey","grey")) 
  z

# HLA-C 
ranges.show <- StringToGRanges(c("chr6-31267604-31268210"))
ranges.show$color <- "darkgrey"
  z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "HLA-C",  features = "HLA-C",
    links = F, 
    extend.upstream = 4000,
    extend.downstream = 0, 
    ymax = 45) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z
  
# GZMA 
ranges.show <- StringToGRanges(c("chr5-55095822-55097308"))
ranges.show$color <- "darkgrey"
    z<- CoveragePlot(
      object = Spike_Specific_CD4, region.highlight	= ranges.show,
      region = "GZMA",  features = "GZMA",
      links = F, 
      extend.upstream = 6339,
      extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#IRF5-- in isg leading edge
ranges.show <- StringToGRanges(c( "chr7-128947044-128948115"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                   region = "chr7-128947044-128948115",  features = "IRF5",
                   links = F, 
                   extend.upstream = 10000,
                   extend.downstream = 4000) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z
  


#BTG1
ranges.show <- StringToGRanges(c("chr12-92184739-92185477", "chr12-92185736-92186423"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(object = Spike_Specific_CD4, region.highlight	= ranges.show,
                   region = "BTG1",
                   links = F, annotation = T,
                   extend.upstream = 0,
                   extend.downstream = 50000, ymax = 25) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#STAT1
ranges.show <- StringToGRanges(c("chr2-191Subject5158-191014895"))
ranges.show$color <- "darkgrey"
DefaultAssay(Spike_Specific_CD4) <- "peaks"
z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "STAT1",  features = "STAT1",
    links = F, annotation = T,
    extend.upstream = 0,
    extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(infect,vax)) & scale_color_manual(values = c("darkgrey","grey")) 
z
#TBX21
ranges.show <- StringToGRanges(c("chr17-47735487-47735832"))
ranges.show$color <- "darkgrey"
z<- CoveragePlot(
    object = Spike_Specific_CD4, region.highlight	= ranges.show,
    region = "TBX21",  features = "TBX21",
    links = F, annotation = T,
    extend.upstream = 0,
    extend.downstream = 0) 
z<- z & scale_fill_manual(values = c(vax,infect)) & scale_color_manual(values = c("darkgrey","grey")) 
z

#--------Heatmap of DAR--------
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
View(Spike_Specific_CD4)

DefaultAssay(Spike_Specific_CD4) <- "peaks"
AIM.reactive.CD4.sce <- Spike_Specific_CD4 %>% as.SingleCellExperiment() #using peak counts !!!
## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$HTO_maxID)
groups <- colData(AIM.reactive.CD4.sce)[, c("HTO_maxID")]
altExpNames(AIM.reactive.CD4.sce)
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("HTO_maxID")],use.assay.type = "logcounts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "-",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject5")  ~ "first",endsWith(subject, "Subject6")  ~ "first", endsWith(subject, "Subject8")  ~ "first", endsWith(subject, "Subject9")  ~ "first",
    endsWith(subject, "Subject1")  ~ "second",endsWith(subject, "Subject24")  ~ "second",endsWith(subject, "Subject17")  ~ "second",
    endsWith(subject, "Subject7") ~ "second", endsWith(subject, "Subject23") ~ "second",  endsWith(subject, "Subject20") ~ "second", endsWith(subject, "Subject12") ~ "second"))
metaData

# Create DESeq2 object        
pseudo_bulk_CD4 <- assay(pseudo_bulk_CD4) 
cluster_counts <- as.data.frame(as.matrix(pseudo_bulk_CD4[, which(colnames(pseudo_bulk_CD4) %in% metaData$subject_timepoint)])) 
#View(cluster_counts) #just making sure metaData and counts have same info
# Check that all of the row names of the metadata are the same 
#and in the same order as the column names of the counts in order to use as input to DESeq2
all(metaData$subject_timepoint == colnames(cluster_counts))         

DefaultAssay(Spike_Specific_CD4) <- "peaks"
View(Spike_Specific_CD4)
ActivatedT.mxt.peaks <- as.data.frame(Spike_Specific_CD4@assays$peaks@data) #36601 genes
#View(ActivatedT.mxt.peaks)
ActivatedT.mxt.peaks <- as.data.frame(ActivatedT.mxt.peaks)

View(top.da.peak)
top.da.peak <- (da_peaks_sig[da_peaks_sig$p_val < 0.005, ])
DAR <- da_peaks_sig$query_region

cluster_counts <- as.data.frame(cluster_counts)
pseudo.bulk.top <- cluster_counts %>% 
  dplyr::filter(rownames(cluster_counts) %in% DAR) %>% 
  as.matrix
pseudo.bulk.top <- scale(pseudo.bulk.top)

metaData <- metaData[,c("subject_timepoint", "timepoint")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  timepoint = c(`PostSecondDose` = vax, `PostInfection` = infect))

z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = F,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA, 
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        show_rownames = F,
                        #labels_row = genes,
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        annotation_legend = F,
                        annotation_names_col = FALSE,
                        color = colorRampPalette(c("#002F70", "white", "#881415"))(50))#color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))#

ggsave(z, filename = "./Images/Images/TEAseq_Heatmap_DAR.pdf", width = 3, height = 4.5)


#---------Motif------------
BiocManager::install("JASPAR2Subject6", force = T)
library(JASPAR2Subject6)
BiocManager::install("TFBSTools", force = T)
library(TFBSTools)
BiocManager::install("motifmatchr")
library(motifmatchr)

DefaultAssay(Spike_Specific_CD4)

pfm <- getMatrixSet(
  x = JASPAR2Subject6,
  opts = list(species = 9606) # 9606 is the species code for human
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(Spike_Specific_CD4),
  pwm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)

Motifs(Spike_Specific_CD4) <- motif

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(Spike_Specific_CD4, idents = c("PostInfection", "PostVaccination"), assay = "peaks")

top.da.peak <- (da_peaks_sig[da_peaks_sig$p_val < 0.005, ])
top.da.peak <- top.da.peak$query_region
# match the overall GC content in the peak set
meta.feature <- GetAssayData(Spike_Specific_CD4, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)

vax_peaks <- subset(vax_open, p_val < 0.005)
enriched.motifs_vax <- FindMotifs(
  object = Spike_Specific_CD4,
  features = vax_peaks$query_region, 
  background=peaks.matched
)

MotifPlot(
  object = Spike_Specific_CD4,
  motifs = head(rownames(enriched.motifs_vax))
)
#View(enriched.motifs_vax)
enriched.motifs_vax$diffexpressed <- "NO"
enriched.motifs_vax$diffexpressed[enriched.motifs_vax$p.adjust < 0.05] <- "VAX"

vax_motifs <- subset(enriched.motifs_vax, subset = diffexpressed == "VAX" )
top <- vax_motifs$motif.name
enriched.motifs_vax_plot <- enriched.motifs_vax
rownames(enriched.motifs_vax_plot) <- NULL
enriched.motifs_vax_plot <- column_to_rownames(enriched.motifs_vax_plot, var = "motif.name" )
top
top <- c("EGR4","HES1","HES2","HEY1","KLF2","SP1","KLF5","KLF4")
plot2<- ggplot(data=enriched.motifs_vax_plot,aes(x=fold.enrichment, y=-log10(p.adjust), col = diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ scale_y_continuous(limits = c(0,8)) +scale_x_continuous(limits = c(0,3))+ 
  scale_color_manual(values = c("lightgrey",vax))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "fold enrichment over background")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_Motif_vax.pdf", width = 4, height =6)

#read in output from ChEA3
df<-read.delim("~/S3data/Gene_sets/vax_imprint_predictedMotif_ChEA3.tsv",sep="\t")
#View(df)
ChEA3.output <- df[,c("TF","Score")]
motifs_vax <- subset(enriched.motifs_vax, diffexpressed == "VAX")
#View(enriched.motifs_vax)
df <- merge(motifs_vax, ChEA3.output, by.x = "motif.name", by.y = "TF")
#View(df)
df <- column_to_rownames(df, var = "motif.name" )
plot2 <- ggplot(df, aes(x = Score, y = -log10(p.adjust), col = diffexpressed)) + geom_point()+theme_classic()+ 
   scale_color_manual(values = c(vax))
top <- c("NRF1" ,"EGR1" ,"EGR4","HES1","HES2","HEY1","KLF2","SP1","KLF5","KLF4", 
         "KLF16","ZBTB14","KLF11","ZNF148","SP9", "SP8","SP4","SP3", "EGR3",
         "ZNF740")
z <- LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+
  NoLegend()+labs(x = "ChEA3 Score", y = "-log10 adjusted p value")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) +theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_Motif_vax_ChEA3_cor.pdf", width = 6, height =6)

chea_corr <- subset(df, diffexpressed == "VAX")
cor.test(-log10(chea_corr$p.adjust), chea_corr$Score, method = 'pearson')



infect_peaks <- subset(infect_open, p_val < 0.005)
enriched.motifs_infect <- FindMotifs(
  object = Spike_Specific_CD4,
  features = infect_peaks$query_region,
  background=peaks.matched
)
#View(enriched.motifs_infect)
MotifPlot(
  object = Spike_Specific_CD4,
  motifs = head(rownames(enriched.motifs_infect))
)
enriched.motifs_infect$diffexpressed <- "NO"
enriched.motifs_infect$diffexpressed[enriched.motifs_infect$p.adjust < 0.05] <- "INFECT"
infect_motifs <- subset(enriched.motifs_infect, subset = diffexpressed == "INFECT" )
top <- infect_motifs$motif.name
enriched.motifs_infect_plot <- enriched.motifs_infect
rownames(enriched.motifs_infect_plot) <- NULL
enriched.motifs_infect_plot <- column_to_rownames(enriched.motifs_infect_plot, var = "motif.name" )
View(enriched.motifs_infect_plot)

top <- c("IRF8","IRF9","IRF2","IRF7","STAT1::STAT2","BATF3", "PRDM1",          
         "IRF4","IRF1","DBP", "MAF","GRHL1","TEF")
plot2<- ggplot(data=enriched.motifs_infect_plot, aes(x=fold.enrichment, y=-log10(p.adjust), col = diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ scale_y_continuous(limits = c(0,8)) +scale_x_continuous(limits = c(0,3))+ 
  scale_color_manual(values = c(infect,"lightgrey"))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6, ynudge = -0.2, xnudge = -0.2 )+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "fold enrichment over background")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_Motif_infect.pdf", width = 4, height =6)

top <- c("IRF8","IRF9","IRF2","IRF7","STAT1::STAT2","BATF3", "PRDM1",          
         "IRF4","IRF1")
plot2<- ggplot(data=enriched.motifs_infect_plot, aes(x=percent.background, y=percent.observed, col = diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+  scale_y_continuous(limits = c(0.1,100)) +scale_x_continuous(limits = c(0.1,100))+ 
  scale_color_manual(values = c(infect,"lightgrey"))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6, ynudge = 0.4, xnudge = -0.5 )+
  NoLegend()+labs(y = "% of target with motif", x = "% of background with motif")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_Motif_infect_percentObserved.pdf", width = 5, height =5)


#read in output from ChEA3
df<-read.delim("~/S3data/Gene_sets/infect_imprint_predictedMotif_ChEA3.tsv",sep="\t")
#View(df)
ChEA3.output <- df[,c("TF","Score")]
#View(enriched.motifs_infect)
motifs_infect <- subset(enriched.motifs_infect, diffexpressed == "INFECT")
df <- merge(motifs_infect, ChEA3.output, by.x = "motif.name", by.y = "TF")
#View(df)
df <- column_to_rownames(df, var = "motif.name" )
plot2 <- ggplot(df, aes(x = Score, y = -log10(p.adjust), col = diffexpressed)) + geom_point()+theme_classic()+ 
  scale_color_manual(values = c(infect,"lightgrey"))
top <- c("IRF8","IRF9","IRF2","IRF7",          
         "IRF4","IRF1","IRF5","IRF3")
z <- LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+
  NoLegend()+labs(x = "ChEA3 Score", y = "-log10 adjusted p value")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) +theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))

z
ggsave(z, filename = "~/S3data/Images/Images/TEAseq_Motif_infect_ChEA3_cor.pdf", width = 6, height =6)

chea_corr <- subset(df, diffexpressed == "INFECT")
cor.test(-log10(chea_corr$p.adjust), chea_corr$Score, method = 'pearson')
