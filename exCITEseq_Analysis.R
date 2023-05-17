library(pheatmap)
#BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
#install.packages('Matrix.utils')
#library(Matrix.utils)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("gplots")
library(gplots)
library(colorspace)
#install.packages("tidyverse")
library(tidyverse)
#BiocManager::install("limma")
library(limma)
#install.packages("dplyr")
library(dplyr)
#devtools::install_github("ncborcherding/scRepertoire@dev", force = TRUE)
library(scRepertoire)
library(readxl)
library(DESeq2)
library(SummarizedExperiment)
#if (!require("BiocManager", quietly = TRUE))
#BiocManager::install("GSEABase")
library(GSEABase)
#BiocManager::install("GSVA")
library(GSVA)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(Seurat)
library(RColorBrewer)
library(viridis)
#BiocManager::install("scuttle")
library(scuttle)
#BiocManager::install("scater")
library(scater)
#BiocManager::install("limma")
library(limma)
install.packages("SCpubr")
library(SCpubr)

sessionInfo()
setwd("~/S3data")

#-------colors------
pre_vax <- "#fdc472"
post_vax <- "#e59628"
pre_infect <- "#b5b2f1"
post_infect <- "#6667a9"
vax <- "#ffc471"
infect <- "#b5b2f1"
breakthru <- "#D96F6F"
colorblind_vector <- colorRampPalette((c("#0D0887FF", "#7301A8FF", 
                                         "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                        "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
                                                  
                                                    
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

#-------------------------Read in AIM.integrated.T.clean object--------------

AIM.integrated.T.clean.CD4 <- readRDS(file = "~/S3data/Saved_Objects/20221214_AIM.integrated.CD4T.rds")
#jumping back in
#AIM.reactive.CD4 <- readRDS(file = "~/S3data/Saved_Objects/20221214_AIM.reactive.CD4.RDS")

#-----------------------Split CD4 UMAP by condition-----------------------
#Split UMAP by condition -- stim and unstim
DimPlot(AIM.integrated.T.clean.CD4, group.by = "seurat_clusters", label=TRUE, split.by="sample")

Idents(AIM.integrated.T.clean.CD4) <- "seurat_clusters"
ClusterBySample <- table(Idents(AIM.integrated.T.clean.CD4), AIM.integrated.T.clean.CD4$orig.ident) %>% as.data.frame() 
ggplot(ClusterBySample, aes(fill=Var2, x=Var1, y = Freq)) + 
  geom_bar(stat = "identity")+theme_classic()+labs(y = "Count", x = "Cluster")
#DOWN SAMPLE and Split UMAP by stim and unstim
Idents(AIM.integrated.T.clean.CD4) <-"integrated_cluster.IDs"
p2 <-DimPlot(AIM.integrated.T.clean.CD4,label=TRUE, label.size = 6, raster= FALSE) + NoLegend()+labs(title="")
p2
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.CD4T_UMAP.integrated_cluster.IDs.pdf", width = 14, height = 10)
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.CD4T_UMAP.integrated_cluster.IDs.jpeg", width = 14, height = 10)

#DOWN SAMPLE and Split UMAP by condition -- stim and unstim
AIM.integrated.T.clean.CD4$sample <- factor(AIM.integrated.T.clean.CD4$sample, levels= c("unstim", "stim"))
DefaultAssay(AIM.integrated.T.clean.CD4) <- 'RNA'
AIM.integrated.T.clean.CD4@meta.data <- AIM.integrated.T.clean.CD4@meta.data %>%
  unite('sample_condition', HTO_maxID,sample, remove=FALSE)
Idents(AIM.integrated.T.clean.CD4) <- "sample_condition"
table(Idents(AIM.integrated.T.clean.CD4))
AIM.integrated.T.clean.CD4.sample <- subset(x = AIM.integrated.T.clean.CD4, downsample = 1000)
table(Idents(AIM.integrated.T.clean.CD4.sample))     #corectly downsampled so 49763 
DimPlot(AIM.integrated.T.clean.CD4.sample, group.by = "seurat_clusters", label=FALSE, split.by="sample")+NoLegend()+labs(title="")
x <- DimPlot(AIM.integrated.T.clean.CD4.sample, group.by = "seurat_clusters", label=FALSE, split.by="sample")+NoLegend()+labs(title="")
z <- DimPlot(AIM.integrated.T.clean.CD4.sample, group.by = "integrated_cluster.IDs", label=F, split.by="sample",
             label.box = F, label.color= "black", raster = F)+NoLegend()+labs(title="")
z
#saving AIM T cell UMAP SPLIT BY condition
ggsave(z, filename = "~/S3data/Images/Images/AIMintegrated.CD4T_Unstim_vs_Stim_UMAP.pdf", width = 15, height = 7)
ggsave(z, filename = "~/S3data/Images/Images/AIMintegrated.CD4T_Unstim_vs_Stim_UMAP.jpeg", width = 15, height = 7)

#saving AIM T cell UMAP
Idents(AIM.integrated.T.clean.CD4) <-"integrated_cluster.IDs"
p2 <-DimPlot(AIM.integrated.T.clean.CD4,label=TRUE, label.size = 6, raster= FALSE) + NoLegend()+labs(title="")
p2
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.T_UMAP.integrated_cluster.IDs.pdf", width = 14, height = 10)
ggsave(p2, filename = "~/S3data/Images/Images/AIMintegrated.T_UMAP.integrated_cluster.IDs.jpeg", width = 14, height = 10)

p2 <- DimPlot(AIM.integrated.T.clean.CD4, group.by = "integrated_cluster.IDs", label=FALSE, 
              cols = c( "lightgray", "lightgray","#00A9FF",
                                   "lightgray","#00A9FF","lightgray","lightgray"), raster=FALSE)+labs(title="")+ NoLegend() 
p2 
ggsave(p2, filename = "./Images/Images/AIMintegrated.T_UMAP.AIMreactiveClusters_highlighted.pdf", width = 14, height = 10)
ggsave(p2, filename = "./Images/Images/AIMintegrated.T_UMAP.AIMreactiveClusters_highlighted.jpeg", width = 14, height = 10)
                                   

#-------------------------differential expression analysis stim versus unstim samples --------------
Idents(AIM.integrated.T.clean.CD4)<- "sample"
AIM.integrated.T.clean.sample <- subset(x = AIM.integrated.T.clean.CD4, downsample = 10000)
table(Idents(AIM.integrated.T.clean.sample))     #correctly downsampled so 5000 (~number of unstim in final dataset and using as limiting factor) cells in each condition
DefaultAssay(AIM.integrated.T.clean.sample) <- "RNA"

unstim.stim.DE <- FindMarkers(AIM.integrated.T.clean.sample, ident.1 = "stim", ident.2 = "unstim", test.use="DESeq2", verbose=TRUE) 
unstim.stim.DE$diffexpressed <- "NO"
unstim.stim.DE$diffexpressed[unstim.stim.DE$avg_log2FC > 0.6 & 
                               unstim.stim.DE$p_val_adj < 0.01] <- "UP"
unstim.stim.DE$diffexpressed[unstim.stim.DE$avg_log2FC < -0.6 & 
                               unstim.stim.DE$p_val_adj < 0.01] <- "DOWN"

stim <- subset(unstim.stim.DE, subset = diffexpressed == "UP" )
unstim <- subset(unstim.stim.DE, subset = diffexpressed == "DOWN" ) 
Stim_UP_Unstim_DOWN <- unstim.stim.DE

top <- rbind(stim,unstim)%>% row.names()
plot2<- ggplot(data=unstim.stim.DE, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  geom_hline(yintercept = 2, linetype = "longdash", color = "lightgrey")+
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+
  scale_color_manual(values = c("grey","darkred"))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 5)+labs(title="Stimulated versus Unstimulated")+
  NoLegend()+xlim(limits=c(-2,2))+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme()

z 
ggsave(z, filename = "./Images/Images/Unstim_Stim_Volcano.pdf", width = 8, height =7)
ggsave(z, filename = "./Images/Images/Unstim_Stim_Volcano.jpeg", width = 8, height =7)
write.csv(Stim_UP_Unstim_DOWN,"./Gene_sets/Stim_UP_Unstim_DOWN.csv", row.names = TRUE)
write.csv(stim,"./Gene_sets/Stim_UP_versus_Unstim.csv", row.names = TRUE)
#---------Metascape analysis of genes enriched in stim versus unstim------------------
metasacpe_result <- read_excel("./Metascape_analyses/METAscape_Unstim_versus_Stim/metascape_result.xlsx", 
                               sheet = "Enrichment")
View(metasacpe_result)
subset <- metasacpe_result %>%
  filter(LogP <= -16.24, str_ends(GroupID, "Summary"))

z<- ggplot(subset, aes(x=-LogP, y=reorder(Description, -LogP), fill = -LogP))+geom_col()+ 
  scale_fill_gradient(low = "lightgrey",high = "darkred",limits=c(0,60))+
  theme(legend.key = element_rect(fill=NA)) + 
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1),
        axis.text=element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15))+labs(y="",x="")
z
ggsave(z, filename = "./Images/Images/EnrichedTerms_Stimulated.jpeg", width = 10, height =8)
ggsave(z, filename = "./Images/Images/EnrichedTerms_Stimulated.pdf", width = 10,height =8)
write.csv(subset,"~/S3data/Export_csv/METAscape_Stim_enriched_terms.csv", row.names = TRUE)

#---------find differentially expressed genes between AIM reactive and non reactive cells------------------
AIM_Reactive_and_Unreactive_CD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "CD4+ T" |
                                            integrated_cluster.IDs == "AIM Reactive CD4+ T")
AIM_Reactive_and_Unreactive_CD4_pre_boost <- subset(AIM_Reactive_and_Unreactive_CD4, timepoint == "7moPostVax")
Idents(AIM_Reactive_and_Unreactive_CD4_pre_boost)<- "integrated_cluster.IDs"
table(AIM_Reactive_and_Unreactive_CD4_pre_boost$integrated_cluster.IDs)
AIM_Reactive_and_Unreactive_CD4_pre_boost.sample <- subset(x = AIM_Reactive_and_Unreactive_CD4_pre_boost, downsample = 1000)
table(AIM_Reactive_and_Unreactive_CD4_pre_boost.sample$integrated_cluster.IDs)
DefaultAssay(AIM_Reactive_and_Unreactive_CD4_pre_boost.sample) <- "RNA"

AIM.reactivity.CD4 <- FindMarkers(AIM_Reactive_and_Unreactive_CD4_pre_boost.sample, ident.1 = "AIM Reactive CD4+ T", ident.2 = "CD4+ T", test.use="DESeq2") 
AIM.reactivity.CD4$diffexpressed <- "NO"
AIM.reactivity.CD4$diffexpressed[AIM.reactivity.CD4$avg_log2FC > 1 & 
                                   AIM.reactivity.CD4$p_val_adj < 0.01] <- "UP"
AIM.reactivity.CD4$diffexpressed[AIM.reactivity.CD4$avg_log2FC < 0 & 
                                   AIM.reactivity.CD4$p_val_adj < 0.01] <- "DOWN"

AIM.CD4 <- subset(AIM.reactivity.CD4, subset = diffexpressed == "UP" )
CD4 <- subset(AIM.reactivity.CD4, subset = diffexpressed == "DOWN" ) 
AIMreactive_UP_NONreactive_DOWN_CD4 <- AIM.reactivity.CD4

top <- AIM.CD4 %>% row.names() 
View(top)
plot2<- ggplot(data=AIM.reactivity.CD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 2, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-2, 2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 2)+ theme_classic()+
  scale_color_manual(values = c("black","grey","#00A9FF"))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 6, max.overlaps = 18)+
  labs(title="CD4+ T unreactive versus AIM reactive")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 

z 
ggsave(z, filename = "./Images/Images/CD4_7moPostVax_AIMreactivity_Volcano.pdf", width = 8, height =8)
ggsave(z, filename = "./Images/Images/CD4_7moPostVax_AIMreactivity_Volcano.jpeg", width = 8, height =8)
write.csv(AIMreactive_UP_NONreactive_DOWN_CD4,"./Gene_sets/AIMreactive_UP_NONreactive_DOWN_CD4.csv", row.names = TRUE)
write.csv(AIM.CD4,"./Gene_sets/AIMreactive_UP_vs_NONreactive_CD4.csv", row.names = TRUE)
AIM.reactivity.CD4 <- read.csv(file = "./Gene_sets/AIMreactive_UP_NONreactive_DOWN_CD4.csv", row.names = 1)
View(AIM.reactivity.CD4)

#---------Metascape analysis of differentially expressed genes between AIM reactive and non reactive cells------------------
metasacpe_result <- read_excel("./Metascape_analyses/METAscape_AIMreactive_UP_vs_NONreactive_CD4/metascape_result.xlsx", 
                               sheet = "Enrichment")
subset <- metasacpe_result %>%
  filter(LogP <= -20, str_ends(GroupID, "Member"))
View(subset)
subset <- subset %>%
  filter(Description %in% c("Metabolism of RNA",
                            "Coronavirus disease - COVID-19",
                            "ATP metabolic process",
                            "TCR signaling",
                            "Downstream TCR signaling", "SARS-CoV-2 Infection", "Signaling by Interleukins"))

z<- ggplot(subset, aes(x=-LogP, y=reorder(Description, -LogP), fill = -LogP))+geom_col()+ scale_fill_gradient(low = "lightgrey", high = "darkred",limits=c(0,100))+
  theme(legend.key = element_rect(fill=NA)) + 
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1),
        axis.text=element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15))+labs(x="",y="")
z
ggsave(z, filename = "./Images/Images/EnrichedTerms_AIM_reactive_CD4_T.jpeg", width = 10, height =4)
ggsave(z, filename = "./Images/Images//EnrichedTerms_AIM_reactive_CD4_T.pdf", width = 12,height =8)
write.csv(subset,"~/S3data/Export_csv/METAscape_EnrichedTerms_AIM_reactive_CD4_T.csv", row.names = TRUE)



#-----------------Cytokines within AIM reactive cluster----------------------------- 
DefaultAssay(AIM.integrated.T.clean.CD4) <-"RNA"
Cytokines <- c("IFNG","TNF","IL2","IL21","IL10", "IL5", "IL13",  "IL4", "IL17A", "IL6")
DotPlot(AIM.integrated.T.clean.CD4, features = Cytokines, group.by = "integrated_cluster.IDs",  
        scale  = FALSE, cols = c("lightgray", "red4")) +
 labs(x="Cytokines") + RotatedAxis()

DefaultAssay(AIM.integrated.T.clean.CD4) <-"RNA"
AIM.reactive.CD4.pre.boost <- subset(AIM.integrated.T.clean.CD4, timepoint == "7moPostVax")
DotPlot(AIM.reactive.CD4.pre.boost, features = Cytokines, scale  = FALSE, cols = c("lightgray", "red4"), group.by = "PriorCOVID", assay = "RNA") +
  coord_flip() +labs(x="Cytokines") + RotatedAxis() + scale_size(limits = c(0,100))
#rm(AIM.reactive.CD4.pre.boost)
#cytokine dot plot with only AIM reactive CD4s compared to unreactive CD4s
AIM_Reactive_and_Unreactive_CD4_pre_boost <- subset(AIM.integrated.T.clean.CD4, timepoint == "7moPostVax")
DefaultAssay(AIM_Reactive_and_Unreactive_CD4_pre_boost) <-"RNA"
Cytokines <- c("TGFB1","TNF","IFNG","IL2","IL21","IL13","IL10", "IL4","IL5",  "IL17A")
z<-DotPlot(AIM_Reactive_and_Unreactive_CD4_pre_boost, features = Cytokines,  group.by = "integrated_cluster.IDs",
           scale  = FALSE, cols = c("gray", "#6F0000"), scale.by = "size", scale.min = 0) +
  labs(x="Cytokines") + coord_flip() + RotatedAxis() +
  theme(panel.background=element_rect(fill="white"))
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_CD4sONLY_IL_DotPlot.pdf", width = 4, height = 7)

DefaultAssay(AIM_Reactive_and_Unreactive_CD4_pre_boost) <-"RNA"
Cytokines <- c( "GZMB","TNF","IL2","IFNG")
Idents(AIM_Reactive_and_Unreactive_CD4_pre_boost) <- "integrated_cluster.IDs"
z<-DotPlot(AIM_Reactive_and_Unreactive_CD4_pre_boost, features = Cytokines,  group.by  = "integrated_cluster.IDs",
           scale  = FALSE, cols = c("gray", "#6F0000"), scale.by = "size", scale.min = 0) +
  labs(x="Cytokines") +
  theme(panel.background=element_rect(fill="white")) + labs(x="",y="") + coord_flip() +RotatedAxis()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_CD4sONLY_4_cytokine_DotPlot.pdf", width = 4, height = 4.5)
rm(AIM_Reactive_and_Unreactive_CD4_pre_boost)

#-----------------Cytokine analsis deep dive (FIGURE 2)----------------------------- 
#cytokine analsis deep dive
AIM.reactive.CD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "AIM Reactive CD4+ T")
DefaultAssay(AIM.reactive.CD4) <-"RNA"
Cytokines <- c("IFNG","TNF","IL2","IL12A","CXCR3","CCR5","STAT4","TBX21","RUNX3",           #Th1
               "IL4", "IL5", "IL13", "CXCR4","GATA3", "STAT6", "CCR4",              #Th2
               "IL21","IL17A", "RORC","STAT3", "CCR6",                        #Th17
               "IL10" , "IL2RA", "FOXP3", "CCR7","TGFB1",                   #Treg
               "CXCR5"
               )

#FeaturePlot(AIM.reactive.CD4, features = Cytokines)
AIM.reactive.CD4 <-AIM.reactive.CD4 %>% NormalizeData() %>%  ScaleData()
AIM.reactive.CD4 <- RunPCA(AIM.reactive.CD4, approx = FALSE, features = Cytokines) 
DimPlot(AIM.reactive.CD4, reduction = "pca")
ElbowPlot(AIM.reactive.CD4)
AIM.reactive.CD4<- RunUMAP(AIM.reactive.CD4,
        dims = 1:15,
        reduction = 'pca',
        assay = 'RNA',
        reduction.name = 'cytokine_umap')
AIM.reactive.CD4 <- FindNeighbors(AIM.reactive.CD4, dims = 1:15)
AIM.reactive.CD4 <- FindClusters(AIM.reactive.CD4, resolution = 0.4, verbose = FALSE)
#Idents(AIM.reactive.CD4) <- "RNA_snn_res.0.4"
z<- DimPlot(AIM.reactive.CD4,  label = TRUE, label.color = "black",  label.size = 5, reduction = "cytokine_umap", group.by = "RNA_snn_res.0.4") +
  scale_color_brewer(palette = "Paired")  + NoAxes()
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP.pdf", width = 4.5, height = 3.5)

z<- DimPlot(AIM.reactive.CD4,  label = TRUE, label.color = "black",  label.size = 5, split.by = "timepoint") +
  scale_color_brewer(palette = "Paired")  + NoAxes()
z

FeatureScatter(AIM.reactive.CD4, feature1 = "GZMB", feature2 = "PRF1")
subset <- subset(AIM.reactive.CD4, timepoint != "PostBreakThru")
z <- DimPlot(subset, split.by = "highlight", label = T)+NoLegend()+NoAxes()+scale_color_brewer(palette = "Paired")  
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_splitBy_PriorCOVID_noBreakThru.pdf", width = 7, height = 3.5)

dat <- data.frame(subset@meta.data)
dat$UMAP1 <- Embeddings(subset, "cytokine_umap")[,1]
dat$UMAP2 <- Embeddings(subset, "cytokine_umap")[,2]
dat_bg <- dat[,-(which(colnames(dat)=="PriorCOVID"))]
density_plot <- ggplot(dat, aes(x=UMAP1, y=UMAP2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + #ndensity calculates the normalized density for each sample--otherwise density would be affected by the number of cells for each sample, which is variable
  geom_point(data=dat_bg, shape=16, size=0.1, alpha=0.2, color="white") +
  scale_fill_gradientn(colours=viridis::magma(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~PriorCOVID, ncol=2) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))
density_plot
ggsave(density_plot, filename = "./Images/Images/Cytokines_UMAP_splitBy_PriorCOVID_noBreakThru_Density.pdf", width = 7, height = 3.5)
ggsave(density_plot, filename = "./Images/Images/Cytokines_UMAP_splitBy_PriorCOVID_noBreakThru_Density.jpeg", width = 7, height = 3.5)

plot.density2d_level <- dat %>%
  ggplot(aes(x=UMAP1, y=UMAP2)) + 
  stat_density_2d(aes(fill = stat(piece)),
                  geom = "polygon", 
                  n = 100, 
                  bins = 10, 
                  contour = T) + 
  scale_x_continuous(limits=c(-6,11)) +
  scale_y_continuous(limits=c(-6,8)) +
  facet_wrap(~PriorCOVID, ncol=2) + 
  scale_fill_viridis(option = "A", na.value = "transparent") + theme_classic()
plot.density2d_level
ggsave(plot.density2d_level, filename = "./Images/Images/Cytokines_UMAP_splitBy_PriorCOVID_noBreakThru_Density.pdf", width = 7, height = 3.5)
ggsave(plot.density2d_level, filename = "./Images/Images/Cytokines_UMAP_splitBy_PriorCOVID_noBreakThru_Density.jpeg", width = 7, height = 3.5)

#looking at expression of chemokines and cytokines 
z <- FeaturePlot(AIM.reactive.CD4, features = c("IFNG","TNF","IL4","IL13","IL21","IL10"),ncol = 6,
            cols = c("grey","darkred"), combine = TRUE) & NoLegend() & NoAxes()
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_FeaturePlot1.pdf", width = 10.5, height = 2.5)

z <- FeaturePlot(AIM.reactive.CD4, features = c("IFNG","TNF","IL4","IL13","IL21","IL10"),ncol = 6,
                 cols = c("grey","darkred"), combine = TRUE) & NoAxes()
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_FeaturePlot1_withLegend.pdf", width = 10.5, height = 2.5)

z <- FeaturePlot(AIM.reactive.CD4, features = c("CXCR3","CCR4","CCR6","RORC","FOXP3","CXCR5"),ncol = 6,
            cols = c("grey","darkred"), combine = TRUE) & NoLegend() & NoAxes()
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_FeaturePlot2.pdf", width = 10.5, height = 2.5)

z <- FeaturePlot(AIM.reactive.CD4, features = c("CXCR3","CCR4","CCR6","RORC","FOXP3","CXCR5"),ncol = 6, 
                 cols = c("grey","darkred"), combine = TRUE) & NoAxes()
z
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_FeaturePlot2__withLegend.pdf", width = 10.5, height = 2.5)



DotPlot(AIM.reactive.CD4, features =  c("TCF7","SELL","CCR7","CXCR3","TBX21","GZMB","CCL5",
                                      "CXCR4","CCR6","FOXP3","RORC",
                                      "CXCR5","TNF","IL21","IL2","IFNG","IL13","IL4","IL5","CCR5",
                                      "CTLA4","LAG3","IL10"
                                      ), cols = c("lightgrey","darkred"), assay = "RNA") + coord_flip()

color = brewer.pal(name = "Paired", n = 11) 
all <- prop.table(table(AIM.reactive.CD4$seurat_clusters)) %>% as.data.frame()
all$Freq <- all$Freq*100
colnames(all) <- c("Cluster","% AIM Reactive CD4+ T")
z <- ggplot(all, aes(x = "AIM CD4+ T", y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col(width = 0.5, colour = "white")  + theme_classic() +scale_fill_manual(values = color)
z

AIM.reactive.CD4@meta.data <- AIM.reactive.CD4@meta.data %>%
  unite('subject_timepoint', subject,timepoint, remove=FALSE)
AIM.reactive.CD4@meta.data <- AIM.reactive.CD4@meta.data %>%
  unite('PriorCOVID_timepoint', PriorCOVID,timepoint, remove=FALSE)
AIM.reactive.CD4$PriorCOVID_timepoint <- factor(AIM.reactive.CD4$PriorCOVID_timepoint, 
                                                          levels = c("no_7moPostVax" ,
                                                                      "yes_7moPostVax",
                                                                     "no_1moPostBoost",
                                                                     "yes_1moPostBoost",
                                                                     "no_PostBreakThru"))
AIM.reactive.CD4$timepoint <- factor(AIM.reactive.CD4$timepoint, 
                                                          levels = c("7moPostVax" ,
                                                                     "1moPostBoost",
                                                                      "PostBreakThru"))
ggplot(AIM.reactive.CD4@meta.data, aes(fill=factor(seurat_clusters), x=timepoint,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") + scale_fill_manual(values = color)

vaccine <- subset(AIM.reactive.CD4, PriorCOVID == "no" & timepoint != "PostBreakThru")
vaccine$timepoint <- factor(vaccine$timepoint, levels = c("7moPostVax","1moPostBoost"))
vax <- prop.table(table(vaccine$seurat_clusters, vaccine$timepoint), margin = 2)
vax <- vax[,1:2] %>% as.data.frame()
vax$Freq <- vax$Freq*100
colnames(vax) <- c("Cluster","Time Point", "% AIM Reactive CD4+ T")
z <- ggplot(vax, aes(x=`Time Point`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
   geom_col( width = 0.5) +
  geom_area(aes(x = c("7moPostVax" = 1.25, "1moPostBoost" = 1.75) [`Time Point`]), 
            colour = "white", alpha = 0.2,
             outline.type = "both") + 
  theme_classic() + scale_fill_manual(values = color) +
  theme(panel.background=element_rect(fill="white"), 
       axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
       axis.title.y.left = element_text(size = 20))+labs(x="")
z
ggsave(z, filename = "./Images/Images/VaccineImprint_pre_post_ClusterAbundance.pdf", width = 4, height = 4)
rm(vaccine)

infection <- subset(AIM.reactive.CD4, PriorCOVID == "yes" )
infect <- prop.table(table(infection$seurat_clusters, infection$timepoint), margin = 2)
infect <- infect[,1:2] %>% as.data.frame()
infect$Freq <- infect$Freq*100
colnames(infect) <- c("Cluster","Time Point", "% AIM Reactive CD4+ T")
z <- ggplot(infect, aes(x=`Time Point`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col( width = 0.5) +
  geom_area(aes(x = c("7moPostVax" = 1.25, "1moPostBoost" = 1.75) [`Time Point`]), 
            colour = "white", alpha = 0.2,
            outline.type = "both") + 
  theme_classic() + scale_fill_manual(values = color) +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="")
z
ggsave(z, filename = "./Images/Images/InfectionImprint_pre_post_ClusterAbundance.pdf", width = 4, height = 4)
rm(infection)

pre <- subset(AIM.reactive.CD4, timepoint == "7moPostVax" )
pre <- prop.table(table(pre$seurat_clusters, pre$PriorCOVID), margin = 2)%>% as.data.frame()
pre$Freq <- pre$Freq*100
colnames(pre) <- c("Cluster","PriorCOVID", "% AIM Reactive CD4+ T")
z <- ggplot(pre, aes(x=`PriorCOVID`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col(width = 0.75, colour = "white") +
  #geom_area(aes(x = c("no" = 1.25, "yes" = 1.75)[`PriorCOVID`]), 
           # position = "fill", colour = "white", alpha = 0.4,
            #outline.type = "both") + 
  theme_classic() + NoLegend() + scale_fill_brewer(palette = "Paired") +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="Prior COVID?", title = "Pre booster")
z
ggsave(z, filename = "./Images/Images/PreBoost_Infection_Vaccine_ClusterAbundance.pdf", width = 3, height = 4)
rm(pre)

post <- subset(AIM.reactive.CD4, timepoint == "1moPostBoost" )
post <- prop.table(table(post$seurat_clusters, post$PriorCOVID), margin = 2)%>% as.data.frame()
post$Freq <- post$Freq*100
colnames(post) <- c("Cluster","PriorCOVID", "% AIM Reactive CD4+ T")
z <- ggplot(post, aes(x=`PriorCOVID`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col(width = 0.75, colour = "white") +
  #geom_area(aes(x = c("no" = 1.25, "yes" = 1.75)[`PriorCOVID`]), 
  # position = "fill", colour = "white", alpha = 0.4,
  #outline.type = "both") + 
  theme_classic() + NoLegend() + scale_fill_brewer(palette = "Paired") +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="Prior COVID?", title = "Post booster")
z
ggsave(z, filename = "./Images/Images/PostBoost_Infection_Vaccine_ClusterAbundance.pdf", width = 3, height = 4)
rm(post)

breakthru <- subset(AIM.reactive.CD4, PriorCOVID == "no" & timepoint != "7moPostVax")
breakthru$timepoint <- factor(breakthru$timepoint, levels = c("1moPostBoost", "PostBreakThru"))
vax <- prop.table(table(breakthru$seurat_clusters, breakthru$timepoint), margin = 2)
vax <- vax[,1:2] %>% as.data.frame()
vax$Freq <- vax$Freq*100
colnames(vax) <- c("Cluster","Time Point", "% AIM Reactive CD4+ T")
z <- ggplot(vax, aes(x=`Time Point`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col( width = 0.5) +
  geom_area(aes(x = c("1moPostBoost" = 1.25, "PostBreakThru" = 1.75) [`Time Point`]), 
            colour = "white", alpha = 0.2,
            outline.type = "both") + 
  theme_classic() + scale_fill_manual(values = color) +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="")
z
ggsave(z, filename = "./Images/Images/Breakthru_VaccinePrimed_ClusterAbundance.pdf", width = 4, height = 4)
rm(breakthru)

breakthru <- subset(AIM.reactive.CD4, timepoint != "7moPostVax" & PriorCOVID == "yes" |
                      timepoint == "PostBreakThru")
breakthru$timepoint <- factor(breakthru$timepoint, levels = c("1moPostBoost", "PostBreakThru"))
View(breakthru@meta.data)
vax <- prop.table(table(breakthru$seurat_clusters, breakthru$timepoint), margin = 2)
vax <- vax[,1:2] %>% as.data.frame()
vax$Freq <- vax$Freq*100
colnames(vax) <- c("Cluster","Time Point", "% AIM Reactive CD4+ T")
z <- ggplot(vax, aes(x=`Time Point`, y=`% AIM Reactive CD4+ T`, fill = Cluster, group = Cluster)) +
  geom_col( width = 0.5) +
  #geom_area(aes(x = c("1moPostBoost" = 1.25, "PostBreakThru" = 1.75) [`Time Point`]), 
   #         colour = "white", alpha = 0.2,
   #         outline.type = "both") + 
  theme_classic() + scale_fill_manual(values = color) +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="")
z
ggsave(z, filename = "./Images/Images/Breakthru_InfectionPrimed_ClusterAbundance.pdf", width = 4, height = 4)
rm(breakthru)

DefaultAssay(AIM.reactive.CD4) <- "RNA"
pbmc.markers <- FindAllMarkers(AIM.reactive.CD4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(pbmc.markers)
top3 <- pbmc.markers %>% subset(`p_val_adj` <= 0.05)
all <- top3
top3 <- all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
AIM.reactive.CD4.cytokines.avg <- AverageExpression(AIM.reactive.CD4, return.seurat = TRUE)
DefaultAssay(AIM.reactive.CD4.cytokines.avg) <- "RNA"

brbg <- brewer.pal(11, "BrBG")
cols <- c(colorRampPalette(c(brbg[1], brbg[6]))(100), 
          colorRampPalette(c(brbg[6], brbg[11]))(100)[-1])

z<-DoHeatmap(AIM.reactive.CD4.cytokines.avg, features = top3$genes, draw.lines = FALSE,
             group.colors = color, angle = 0, label = T) + #scale_fill_viridis_c(direction = 1, option = "mako")+
  theme(text =  element_text(size = 15)) +   scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = F)
z

ggsave(z, filename = "./Images/Images/Cytokines_HeatMap_top3genes.pdf", width = 5, height = 7)
ggsave(z, filename = "./Images/Images/Cytokines_HeatMap_top3genes.jpeg", width = 5, height = 7)

AIM.reactive.CD4.cytokines.avg <- AverageExpression(AIM.reactive.CD4, return.seurat = TRUE, assays = "ADT", slot = "scale.data")
colors <- brewer.pal(name = "Paired", n = 9)
rownames(AIM.reactive.CD4.cytokines.avg[["ADT"]])
z<-DoHeatmap(AIM.reactive.CD4.cytokines.avg, features = c("CD27.1","CD45RA","CD197","CD127", "CD62L","CXCR5.1",
                                                          "CD183"), draw.lines = FALSE,
             group.colors = colors, angle = 0) + scale_fill_viridis_c(direction = 1, option = "inferno")
z
AIM.reactive.CD4$subject_timepoint <- as.factor(AIM.reactive.CD4$subject_timepoint)
levels(AIM.reactive.CD4.cytokines$subject_timepoint )
AIM.reactive.CD4.cytokines$subject_timepoint <- factor(AIM.reactive.CD4.cytokines$subject_timepoint, levels = c(
  "Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru",
  "Subject2_7moPostVax","Subject2_1moPostBoost" ,"Subject2_PostBreakThru", 
  "Subject3_7moPostVax" ,"Subject3_1moPostBoost" ,
  "Subject4_7moPostVax", "Subject4_1moPostBoost","Subject4_PostBreakThru",
  "Subject5_7moPostVax" ,   "Subject5_1moPostBoost"  , "Subject5_PostBreakThru",
  "Subject6_7moPostVax" , "Subject6_1moPostBoost" ,"Subject6_PostBreakThru", 
  "Subject7_7moPostVax","Subject7_1moPostBoost" ,
  "Subject8_7moPostVax"  ,  "Subject8_1moPostBoost",  
  "Subject9_7moPostVax","Subject9_1moPostBoost",
  "Subject10_7moPostVax"   ,"Subject10_1moPostBoost"  ,
   "Subject11_7moPostVax"   ,"Subject11_1moPostBoost" ,
  "Subject12_7moPostVax"   ,"Subject12_1moPostBoost" ,
  "Subject13_7moPostVax"   ,   "Subject13_1moPostBoost" , 
  "Subject14_7moPostVax" ,    "Subject14_1moPostBoost"
))
order <- levels(AIM.reactive.CD4.cytokines$subject_timepoint)
a <- prop.table(table(AIM.reactive.CD4.cytokines$seurat_clusters, 
                                           AIM.reactive.CD4.cytokines$subject_timepoint), margin = 2) %>% as.data.frame()
write.table(a, file = "~/S3data/Export_csv/CytokineClusterAbundancePerSample.csv",sep=",",row.names = F, col.names = T)

AIM.reactive.CD4.postboost <- subset(AIM.reactive.CD4, timepoint == "1moPostBoost")
Idents(AIM.reactive.CD4.postboost) <- "PriorCOVID"
Cytokines <- c("IL2","TNF","IFNG", "GZMB", "PRF1","GZMH")
z<-DotPlot(AIM.reactive.CD4.postboost, features = Cytokines,  group.by  = "PriorCOVID",
           scale  = F, cols = c("gray", "#6F0000"), scale.by = "size", scale.min = 0) +
  labs(x="Cytokines") +
  theme(panel.background=element_rect(fill="white")) + labs(x="",y="")
z + coord_flip()
VlnPlot(AIM.reactive.CD4.postboost, features = Cytokines,  
        group.by  = "PriorCOVID", log = T )

z<- DimPlot(AIM.reactive.CD4, group.by = "cloneType", raster=F, order = c("Hyperexpanded (0.1 < X <= 1)", 
                                                                          "Large (0.01 < X <= 0.1)" , 
                                                                          "Medium (0.001 < X <= 0.01)"  , 
                                                                          "Small (1e-04 < X <= 0.001)", NA), size = 3) +
  scale_color_manual(values = colorblind_vector(5),na.value=NA) + 
  theme(plot.title = element_blank())
z
unique(AIM.reactive.CD4$cloneType)
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_CloneType.jpeg", width = 8, height = 4)
ggsave(z, filename = "./Images/Images/Cytokines_UMAP_CloneType.pdf", width = 8, height = 4)

#-----------------Polyfunctionality (FIGURE 2)-------------------
c("INFG","IL2","TNF")

#NaiveCD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "Naive CD4+ T")
#FeatureScatter(NaiveCD4, feature1 = "IFNG", feature2 = "GZMB", group.by = "PriorCOVID") + geom_vline(xintercept = 1) +  geom_hline(yintercept = 1)

PolyFunctional <- AIM.reactive.CD4 
DefaultAssay(PolyFunctional) <- 'RNA'
# FeatureScatter(PolyFunctional, feature1 = "IFNG", feature2 = "GZMB", group.by = "PriorCOVID") + geom_vline(xintercept = 0.25) +  geom_hline(yintercept = 0.25)
PolyFunctional.IFNG <- subset(PolyFunctional, subset = IFNG >0.25) 
PolyFunctional.IFNG@meta.data$IFNG <- "+"
PolyFunctional.IFNG.barcodes <- row.names(PolyFunctional.IFNG@meta.data)

PolyFunctional.IL2 <- subset(PolyFunctional, subset = IL2 > 0.25) 
PolyFunctional.IL2@meta.data$IL2 <- "+"
PolyFunctional.IL2.barcodes <- row.names(PolyFunctional.IL2@meta.data)

PolyFunctional.TNF <- subset(PolyFunctional, subset = TNF > 0.25) 
PolyFunctional.TNF@meta.data$TNF <- "+"
PolyFunctional.TNF.barcodes <- row.names(PolyFunctional.TNF@meta.data)

PolyFunctional.GZMB <- subset(PolyFunctional, subset = GZMB > 0.25) 
PolyFunctional.GZMB@meta.data$GZMB <- "+"
PolyFunctional.GZMB.barcodes <- row.names(PolyFunctional.GZMB@meta.data)

PolyFunctional.IL10 <- subset(PolyFunctional, subset = IL10 > 0.25) 
PolyFunctional.IL10@meta.data$IL10 <- "+"
PolyFunctional.IL10.barcodes <- row.names(PolyFunctional.IL10@meta.data)

PolyFunctional.input <- PolyFunctional@meta.data %>%
  add_column(IFNG = 
               row.names(PolyFunctional@meta.data) %in% PolyFunctional.IFNG.barcodes) 

PolyFunctional.input <- PolyFunctional.input %>%
  add_column(IL2 = 
               row.names(PolyFunctional.input) %in% PolyFunctional.IL2.barcodes) 

PolyFunctional.input <- PolyFunctional.input %>%
  add_column(TNF = 
               row.names(PolyFunctional.input) %in% PolyFunctional.TNF.barcodes)

PolyFunctional.input <- PolyFunctional.input %>%
  add_column(GZMB = 
               row.names(PolyFunctional.input) %in% PolyFunctional.GZMB.barcodes) 

PolyFunctional.input <- PolyFunctional.input %>%
  add_column(IL10 = 
               row.names(PolyFunctional.input) %in% PolyFunctional.IL10.barcodes) 

PolyFunctional.input$IL2[PolyFunctional.input$IL2 == "TRUE"] <- "+"      
PolyFunctional.input$IFNG[PolyFunctional.input$IFNG == "TRUE"] <- "+" 
PolyFunctional.input$TNF[PolyFunctional.input$TNF == "TRUE"] <- "+"
PolyFunctional.input$GZMB[PolyFunctional.input$GZMB == "TRUE"] <- "+"      
PolyFunctional.input$IL10[PolyFunctional.input$IL10 == "TRUE"] <- "+" 


PolyFunctional.input$IL2[PolyFunctional.input$IL2 == "FALSE"] <- "-"      
PolyFunctional.input$IFNG[PolyFunctional.input$IFNG == "FALSE"] <- "-" 
PolyFunctional.input$TNF[PolyFunctional.input$TNF == "FALSE"] <- "-" 
PolyFunctional.input$GZMB[PolyFunctional.input$GZMB == "FALSE"] <- "-"      
PolyFunctional.input$IL10[PolyFunctional.input$IL10 == "FALSE"] <- "-" 

View(PolyFunctional.input.ready)
PolyFunctional.input.ready <- data.frame(PolyFunctional.input$subject, PolyFunctional.input$timepoint, PolyFunctional.input$PriorCOVID,  
                                          PolyFunctional.input$IFNG, PolyFunctional.input$IL2, PolyFunctional.input$TNF, PolyFunctional.input$GZMB)

colnames(PolyFunctional.input.ready) <- c("subject","timepoint","PriorCOVID",  "IFNG", "IL2","TNF", "GZMB")

PolyFunctionality <- PolyFunctional.input.ready %>% 
  group_by(subject,timepoint,PriorCOVID,IFNG,IL2,TNF, GZMB) %>%
 dplyr::summarize(count= n()) %>% dplyr::ungroup() %>% as.data.frame()
colnames(PolyFunctionality)
PolyFunctionality <- as.data.frame(PolyFunctionality)

TotalCells <- aggregate(PolyFunctionality$count, by=list(PolyFunctionality$subject,PolyFunctionality$timepoint), FUN=sum)
colnames(TotalCells) <- c("subject","timepoint","TotalCells")

PolyFunctionality <- merge(PolyFunctionality,TotalCells, by = c("subject","timepoint") )

PolyFunctionality$VALUE <- PolyFunctionality$count/PolyFunctionality$TotalCells 

View(PolyFunctionality)
x<-PolyFunctionality[1:16, 4:7]
PolyFunctionality <- unite(PolyFunctionality,subject_timepoint, subject,timepoint,PriorCOVID,remove=TRUE )
PolyFunctionality$subject_timepoint <- as.factor(PolyFunctionality$subject_timepoint)
n <- levels(PolyFunctionality$subject_timepoint)
n

frame <- data.frame(subject_timepoint=rep(n,ea=NROW(x)),x)
PolyFunctionality <- merge(frame,PolyFunctionality,by.x = 1:5, all.x = TRUE)
PolyFunctionality <- separate(PolyFunctionality,subject_timepoint, sep="_", remove = TRUE, 
                                                      into = c("subject","timepoint","PriorCOVID") )
PolyFunctionality[is.na(PolyFunctionality)] <- 0
PolyFunctionality<- as.data.frame(PolyFunctionality)
View(PolyFunctionality)

check <- PolyFunctionality %>%
  group_by(subject,timepoint) %>%
   dplyr::summarize(SUM = sum(VALUE)) 
View(check)
PolyFunctionality <- dplyr::select(PolyFunctionality, -c('count', 'TotalCells'))
write.csv(PolyFunctionality,"~/S3data/Export_csv/PolyFunctional_CD4_withGZMB.csv", row.names = FALSE)

PolyFunctional_CD4_timepoint_priorCOVID <- read_csv("~/S3data/Export_csv/PolyFunctional_CD4.csv")

#-----------------Activated genes within AIM reactive cluster----------------------------- 
DefaultAssay(AIM.integrated.T.clean.CD4) <-"RNA"
z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features = "MIR155HG", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoAxes()+NoLegend()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_MIR155HG_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_MIR155HG_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features = "TNFRSF9", 
               order = T, cols = c("lightgray","darkred"),raster=FALSE)+NoLegend()+NoAxes()+labs(title = "TNFRSF9 (CD137)")
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNFRSF9_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNFRSF9_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features = "TNFRSF4", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()+labs(title = "TNFRSF4 (OX40)")
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNFRSF4_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNFRSF4_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="TNF", 
               order= T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNF_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNF_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="IFNG", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_IFNG_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_IFNG_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="TFRC", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TFRC_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TFRC_UMAP.jpeg", width = 7, height = 5)


z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="RGS16", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_RGS16_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_RGS16_UMAP.jpeg", width = 7, height = 5)


z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="LTA", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_LTA_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_LTA_UMAP.jpeg", width = 7, height = 5)


z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="ZBTB32", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_ZBTB32_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_ZBTB32_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="TNF", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNF_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_TNF_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="IL2", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_IL2_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_IL2_UMAP.jpeg", width = 7, height = 5)

z<-FeaturePlot(object = AIM.integrated.T.clean.CD4, reduction = 'umap', features ="GZMB", 
               order = T, cols = c("lightgray","darkred"), raster=FALSE)+NoLegend()+NoAxes()
z
ggsave(z, filename = "./Images/Images/AIMintegrated.T_GZMB_UMAP.pdf", width = 7, height = 5)
ggsave(z, filename = "./Images/Images/AIMintegrated.T_GZMB_UMAP.jpeg", width = 7, height = 5)

#---------pseudobulk of experienced and naive in Nonnaive CD4 7 mo post vax-------------
NonnaiveCD4.postBoost <- subset(NonnaiveCD4, timepoint == "7moPostVax")
NonnaiveCD4.postBoost@meta.data <- NonnaiveCD4.postBoost@meta.data %>%
  unite('subject_timepoint', subject,timepoint, remove=FALSE)
CD4.sce <- NonnaiveCD4.postBoost %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(CD4.sce$subject_timepoint)
groups <- colData(CD4.sce)[, c("subject_timepoint")]
CD4.sce <- removeAltExps(CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(CD4.sce, ids = colData(CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ PriorCOVID) 

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("7moPostVax_no","7moPostVax_yes")) #"1moPostBoost_no","PostBreakThru_no","1moPostBoost_yes"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_vax,post_infect))#"#ffa527","#56e39d", "#716be4" )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values =c(post_vax,post_infect))#"#ffa527","#56e39d", "#716be4" )) 
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 3 # filter for genes with at least 10 counts in 3/14 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("PriorCOVID","yes","no"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/NonnaiveCD4_Infect_UP_Vax_DOWN_PreBoost_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")]
write.table(geneRank, file = "~/S3data/Gene_sets/NonnaiveCD4_Infect_UP_Vax_DOWN_PreBoost_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F,sep = "\t")

#---------find differentially expressed genes between experienced and naive in CD4  AIM reactive cells 7 mo post vax-------------
#AIM.reactive.CD4 <- subset(AIM.integrated.T.clean, integrated_cluster.IDs == "AIM Reactive CD4+ T" )
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.7moPostVax <- subset(AIM.reactive.CD4, timepoint == "7moPostVax")

Idents(AIM.reactive.CD4.7moPostVax)<- "PriorCOVID"
table(AIM.reactive.CD4.7moPostVax$PriorCOVID)

PriorCOVID.CD4.7moPostVax <- FindMarkers(AIM.reactive.CD4.7moPostVax, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
PriorCOVID.CD4.7moPostVax$diffexpressed <- "NO"
PriorCOVID.CD4.7moPostVax$diffexpressed[PriorCOVID.CD4.7moPostVax$avg_log2FC > 0.05 & 
                                          PriorCOVID.CD4.7moPostVax$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4.7moPostVax$diffexpressed[PriorCOVID.CD4.7moPostVax$avg_log2FC < -0.05 & 
                               PriorCOVID.CD4.7moPostVax$p_val_adj < 0.05] <- "DOWN"

Experienced.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4.7moPostVax, subset = diffexpressed == "UP" )
Naive.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4.7moPostVax, subset = diffexpressed == "DOWN" )
Experienced_UP_Naive_DOWN_CD4.DESeq.7moPostVax <- PriorCOVID.CD4.7moPostVax
PriorCOVID.CD4.7moPostVax$p_val
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)%>% row.names()
#top <- c("HLA-B","HLA-A","IFITM3","IL32",
#         "REL","CCR7","NINJ1", "NUDT4")
plot2<- ggplot(data=PriorCOVID.CD4.7moPostVax, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(pre_vax,"lightgrey",pre_infect))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 6, max.overlaps = 15)+labs(title="Pre Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_7moPostVax_PriorCOVID_Volcano.pdf", width = 7, height =7)
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_7moPostVax_PriorCOVID_Volcano.jpeg", width = 7, height =7)

#write.csv(Experienced_UP_Naive_DOWN_CD4.DESeq.7moPostVax,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PreBoost.DESeq.csv", row.names = TRUE)
PriorCOVID.CD4.7moPostVax <- read.csv(file = "./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PreBoost.DESeq.csv", row.names = 1)
View(PriorCOVID.CD4.7moPostVax)

top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)
z<- ggplot(top, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(pre_vax,pre_infect))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=after_stat(count)), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,670)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_7moPostVax_PriorCOVID_DEGsummary_BarGraph.jpeg", width = 6, height =6)

#violin plots for IFNG, GZMB, IL2, TNF
DefaultAssay(AIM.reactive.CD4.7moPostVax) <- "RNA"
z <- VlnPlot(AIM.reactive.CD4.7moPostVax, features = c("GZMB"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(pre_vax, pre_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_GZMB_pre.pdf", width = 2.5, height = 2)

z <- VlnPlot(AIM.reactive.CD4.7moPostVax, features = c("IL2"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(pre_vax, pre_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_IL2_pre.pdf",width = 2.5, height = 2)

z <- VlnPlot(AIM.reactive.CD4.7moPostVax, features = c("IFNG"),ncol = 1,
               group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(pre_vax, pre_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_IFNG_pre.pdf", width = 2.5, height = 2)

z <- VlnPlot(AIM.reactive.CD4.7moPostVax, features = c("TNF"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(pre_vax, pre_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_TNF_pre.pdf", width = 2.5, height = 2)


Genes <- rownames_to_column(PriorCOVID.CD4.7moPostVax, "genes")
IFNG <- subset(Genes, genes == "IFNG" )
TNF <- subset(Genes, genes == "TNF" )
GZMB <- subset(Genes, genes == "GZMB" )
IL2 <- subset(Genes, genes == "IL2" )
cyto <- rbind(IFNG, TNF, GZMB, IL2)

ggplot(cyto, aes(x = reorder(genes, -avg_log2FC), y = avg_log2FC, fill = p_val)) + geom_bar(stat = "identity") + 
  scale_fill_continuous(type = "viridis") + theme_classic()

rm(AIM.reactive.CD4.7moPostVax)
#---------find differentially expressed genes between experienced and naive in CD4 AIM reactive cells post boost-------------
Idents(AIM.reactive.CD4)<- "PriorCOVID"
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.1moPostBoost <- subset(AIM.reactive.CD4, timepoint == "1moPostBoost")
table(AIM.reactive.CD4.1moPostBoost$PriorCOVID)
PriorCOVID.CD4 <- FindMarkers(AIM.reactive.CD4.1moPostBoost, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"

Experienced.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Naive.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
Experienced_UP_Naive_DOWN_CD4.DESeq <- PriorCOVID.CD4
#View(Experienced.AIMreactive.CD4.DEseq)
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)%>% row.names()
#top <- c("HLA-B","CCL5","GZMH","NKG7","IFITM1","CCL4","IFITM3","CCL3","GZMB","IFI6","LTB","CST7",
#         "CCR7","DDIT4","NFKB2")
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(post_vax,"lightgrey",post_infect))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+labs(title="Post Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_1moPostBoost_PriorCOVID_Volcano.pdf", width = 6, height =6)
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_1moPostBoost_PriorCOVID_Volcano.jpeg", width = 6, height =6)
#write.csv(Experienced_UP_Naive_DOWN_CD4.DESeq,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost.DESeq.csv", row.names = TRUE)
PriorCOVID.CD4 <- read.csv(file = "./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost.DESeq.csv", row.names = 1)
View(Experienced_UP_Naive_DOWN_CD4.DESeq)
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)
z<- ggplot(top, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(post_vax,post_infect))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=..count..), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,670)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_1moPostBoost_PriorCOVID_DEGsummary_BarGraph.jpeg", width = 6, height =6)

#violin plots for IFNG, GZMB, IL2, TNF
DefaultAssay(AIM.reactive.CD4.1moPostBoost) <- "RNA"
z <- VlnPlot(AIM.reactive.CD4.1moPostBoost, features = c("GZMB"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(post_vax, post_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_GZMB_post.pdf", width = 2.5, height = 2)

z <- VlnPlot(AIM.reactive.CD4.1moPostBoost, features = c("IL2"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(post_vax, post_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_IL2_post.pdf",width = 2.5, height = 2)

z <- RidgePlot(AIM.reactive.CD4.1moPostBoost, features = c("IFNG"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(post_vax, post_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_IFNG_post.pdf", width = 2.5, height = 2)

z <- VlnPlot(AIM.reactive.CD4.1moPostBoost, features = c("TNF"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 7) & 
  scale_fill_manual(values = c(post_vax, post_infect))  & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_TNF_post.pdf", width = 2.5, height = 2)


Genes <- rownames_to_column(PriorCOVID.CD4, "genes")
IFNG <- subset(Genes, genes == "IFNG" )
TNF <- subset(Genes, genes == "TNF" )
GZMB <- subset(Genes, genes == "GZMB" )
IL2 <- subset(Genes, genes == "IL2" )
cyto <- rbind(IFNG, TNF, GZMB, IL2)

ggplot(cyto, aes(x = reorder(genes, -avg_log2FC), y = avg_log2FC, fill = p_val)) + geom_bar(stat = "identity") + 
  scale_fill_continuous(type = "viridis") + theme_classic()

Infection_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Vaccine_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
View(Infection_imprint)
rm(AIM.reactive.CD4.1moPostBoost)
#---------find DEG between experienced and naive in CD4 AIM reactive cells post boost WITH CCR7 cutoff-------------
Idents(AIM.reactive.CD4)<- "PriorCOVID"
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.1moPostBoost <- subset(AIM.reactive.CD4, timepoint == "1moPostBoost")
AIM.reactive.CD4.1moPostBoost <- subset(AIM.reactive.CD4.1moPostBoost, CCR7 > 1)

table(AIM.reactive.CD4.1moPostBoost$PriorCOVID)
PriorCOVID.CD4 <- FindMarkers(AIM.reactive.CD4.1moPostBoost, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"

Experienced.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Naive.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
Experienced_UP_Naive_DOWN_CD4.DESeq <- PriorCOVID.CD4
#View(Experienced.AIMreactive.CD4.DEseq)
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)%>% row.names()
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1.5)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(post_vax,"lightgrey",post_infect))
z <-LabelPoints(plot = plot2, points = top, repel = T, size = 6 )+labs(title="Post Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "./Images/Images/CCR7+AIMreactive_CD4_1moPostBoost_PriorCOVID_Volcano.pdf", width = 6, height =6)
ggsave(z, filename = "./Images/Images/CCR7+AIMreactive_CD4_1moPostBoost_PriorCOVID_Volcano.jpeg", width = 6, height =6)
#write.csv(Experienced_UP_Naive_DOWN_CD4.DESeq,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost.DESeq.csv", row.names = TRUE)
View(Experienced_UP_Naive_DOWN_CD4.DESeq)
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)
z<- ggplot(top, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(post_vax,post_infect))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=..count..), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,670)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/CCR7+AIMreactive_CD4_1moPostBoost_PriorCOVID_DEGsummary_BarGraph.jpeg", width = 6, height =6)

list <- list(rownames(Infection_imprint), rownames(Vaccine_imprint),
             rownames(Experienced.AIMreactive.CD4.DEseq), rownames(Naive.AIMreactive.CD4.DEseq))
names(list) = c( "Infection_imprint", "Vaccine_imprint", 
                "CCR7+_infection","CCR7+_vaccine")
# Pairwise venn diagram
z<-ggvenn(list, c("Infection_imprint","CCR7+_infection"),fill_color = c(post_infect, "grey"),auto_scale = T, 
          set_name_size = 0)   
z
ggsave(z, filename = "./Images/Images/CCR7+_AIM_CD4_PostBoost_overlap_Infection_Imprint.jpeg", width = 4, height =4)

z<- ggvenn(list, c("Vaccine_imprint","CCR7+_vaccine"),fill_color = c( post_vax, "grey"), auto_scale = T, 
           set_name_size = 0)
z
ggsave(z, filename = "./Images/Images/CCR7+_AIM_CD4_PostBoost_overlap_Vax_Imprint.jpeg", width = 4, height =4)


#---------find differentially expressed genes between experienced and naive in Nonnaive CD4 7 mo post vax-------------
NonnaiveCD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "CD4+ T" )
DefaultAssay(NonnaiveCD4) <- "RNA"
NonnaiveCD4.postBoost <- subset(NonnaiveCD4, timepoint == "1moPostBoost")

Idents(NonnaiveCD4.postBoost)<- "PriorCOVID"
table(NonnaiveCD4.postBoost$PriorCOVID)

NonnaiveCD4.postBoost <- FindMarkers(NonnaiveCD4.postBoost, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
NonnaiveCD4.postBoost$diffexpressed <- "NO"
NonnaiveCD4.postBoost$diffexpressed[NonnaiveCD4.postBoost$avg_log2FC > 0.05 & 
                                      NonnaiveCD4.postBoost$p_val_adj < 0.05] <- "UP"
NonnaiveCD4.postBoost$diffexpressed[NonnaiveCD4.postBoost$avg_log2FC < -0.05 & 
                                      NonnaiveCD4.postBoost$p_val_adj < 0.05] <- "DOWN"

#write.csv(NonnaiveCD4.postBoost,"./Gene_sets/NonnaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = TRUE)
NonnaiveCD4.postBoost <- read.csv("./Gene_sets/NonnaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = 1)

Experienced.CD4.DEseq <- subset(NonnaiveCD4.postBoost, subset = diffexpressed == "UP" )
Naive.CD4.DEseq <- subset(NonnaiveCD4.postBoost, subset = diffexpressed == "DOWN" )

top <- rbind(Experienced.CD4.DEseq,Naive.CD4.DEseq)%>% row.names()
plot2<- ggplot(data=NonnaiveCD4.postBoost, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,2)) + #scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(pre_vax,"lightgrey",pre_infect))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 6, max.overlaps = 15)+labs(title="Pre Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z

#subsetting on just Naive CD4 
NaiveCD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "Naive CD4+ T" )
DefaultAssay(NaiveCD4) <- "RNA"
NaiveCD4 <- subset(NaiveCD4, timepoint == "1moPostBoost" )

Idents(NaiveCD4)<- "PriorCOVID"
table(NaiveCD4$PriorCOVID)

NaiveCD4 <- FindMarkers(NaiveCD4, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
NaiveCD4$diffexpressed <- "NO"
NaiveCD4$diffexpressed[NaiveCD4$avg_log2FC > 0.05 & 
                         NaiveCD4$p_val_adj < 0.05] <- "UP"
NaiveCD4$diffexpressed[NaiveCD4$avg_log2FC < -0.05 & 
                         NaiveCD4$p_val_adj < 0.05] <- "DOWN"

#write.csv(NaiveCD4,"./Gene_sets/NaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = TRUE)
NaiveCD4 <- read.csv("./Gene_sets/NaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = 1)

NaiveCD4.Infect.DEseq <- subset(NaiveCD4, subset = diffexpressed == "UP" )
NaiveCD4.Vax.DEseq <- subset(NaiveCD4, subset = diffexpressed == "DOWN" )

top <- rbind(NaiveCD4.Infect.DEseq,NaiveCD4.Vax.DEseq)%>% row.names()
plot2<- ggplot(data=NaiveCD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1)+ theme_classic()+ #scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(pre_vax,"lightgrey",pre_infect))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 6, max.overlaps = 15)+labs(title="Post Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z

top <- rbind(Experienced.DEseq,Naive.DEseq)
z<- ggplot(top, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(pre_vax,pre_infect))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=after_stat(count)), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,670)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/NaiveCD4_PostBoost_PriorCOVID_DEGsummary_BarGraph.jpeg", width = 6, height =6)


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

PriorCOVID.CD4 <- read.csv(file = "~/S3data/Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost.DESeq.csv", row.names = 1)
Infection_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
Vaccine_imprint <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )


list <- list(ISG, prolif, 
             rownames(Infection_imprint), rownames(Vaccine_imprint),
             rownames(Experienced.DEseq), rownames(Naive.DEseq))
names(list) = c("ISG","Prolif",
                "Infection_imprint", "Vaccine_imprint", 
                "Naive_infection_gex_up","Naive_vaccine_gex_up")
# Pairwise venn diagram
z<-ggvenn(list, c("Infection_imprint","Naive_infection_gex_up"),fill_color = c(post_infect, "grey"),auto_scale = T, 
          set_name_size = 0)   
z
ggsave(z, filename = "./Images/Images/NaiveCD4_PostBoost_overlap_Infection_Imprint.jpeg", width = 4, height =4)

z<- ggvenn(list, c("Vaccine_imprint","Naive_vaccine_gex_up"),fill_color = c( post_vax, "grey"), auto_scale = T, 
           set_name_size = 0)
z
ggsave(z, filename = "./Images/Images/NaiveCD4_PostBoost_overlap_Vax_Imprint.jpeg", width = 4, height =4)


list <- list(rownames(Infection_imprint), rownames(Vaccine_imprint),
             rownames(Experienced.CD4.DEseq), rownames(Naive.CD4.DEseq))
names(list) = c( "Infection_imprint", "Vaccine_imprint", 
                 "NonnaiveCD4_infect","NonnaiveCD4_vax")

# Pairwise venn diagram
ggvenn(list, c("Infection_imprint","NonnaiveCD4_infect"),fill_color = c( infect, "lightblue"), auto_scale = T)
#auto_scale = T, set_name_size = 0)          
ggvenn(list, c("Vaccine_imprint","NonnaiveCD4_vax"),fill_color = c(vax, "lightblue"), auto_scale = T) 
#  auto_scale = T, set_name_size = 0)

Infection_imprint <-  Infection_imprint%>% 
  add_column(overlap_sign =
               rownames(Infection_imprint) %in% rownames(Infection_imprint))
Infection_imprint$cluster <- "AIM Reactive CD4+ T"

Experienced.CD4.DEseq <-  Experienced.CD4.DEseq%>% 
  add_column(overlap_sign =
               rownames(Experienced.CD4.DEseq) %in% rownames(Infection_imprint))
Experienced.CD4.DEseq$cluster <- "Nonnaive CD4+ T"

NaiveCD4.Infect.DEseq <-  NaiveCD4.Infect.DEseq%>% 
  add_column(overlap_sign =
               rownames(NaiveCD4.Infect.DEseq) %in% rownames(Infection_imprint))
NaiveCD4.Infect.DEseq$cluster <- "Naive CD4+ T"

NonnaiveCD4 <- rbind(Infection_imprint,Experienced.CD4.DEseq, NaiveCD4.Infect.DEseq)
z<- ggplot(NonnaiveCD4, aes(fill=factor(overlap_sign), x=cluster,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Proportion", x = "Cluster") + 
  scale_fill_manual(values = c("lightgray",post_infect))+ NoLegend() + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 20)) 
z
ggsave(z, filename = "./Images/Images/NonnaiveCD4_PostBoost_Infect_Overlap_DEGsummary_BarGraph.jpeg", width = 4, height =4.5)


Vaccine_imprint <-  Vaccine_imprint%>% 
  add_column(overlap_sign =
               rownames(Vaccine_imprint) %in% rownames(Vaccine_imprint))
Vaccine_imprint$cluster <- "AIM Reactive CD4+ T"

Naive.CD4.DEseq <-  Naive.CD4.DEseq%>% 
  add_column(overlap_sign =
               rownames(Naive.CD4.DEseq) %in% rownames(Vaccine_imprint))
Naive.CD4.DEseq$cluster <- "Nonnaive CD4+ T"

NaiveCD4.Vax.DEseq <-  NaiveCD4.Vax.DEseq%>% 
  add_column(overlap_sign =
               rownames(NaiveCD4.Vax.DEseq) %in% rownames(Vaccine_imprint))
NaiveCD4.Vax.DEseq$cluster <- "Naive CD4+ T"


NonnaiveCD4 <- rbind(Vaccine_imprint,Naive.CD4.DEseq,NaiveCD4.Vax.DEseq)
z<- ggplot(NonnaiveCD4, aes(fill=factor(overlap_sign), x=cluster,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Proportion", x = "") + 
  scale_fill_manual(values = c("lightgray",post_vax)) + NoLegend() + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 20)) 
z
ggsave(z, filename = "./Images/Images/NonnaiveCD4_PostBoost_Vax_Overlap_DEGsummary_BarGraph.jpeg", width = 4, height =4.5)


#subsetting on just unstim
NonnaiveCD4 <- subset(AIM.integrated.T.clean.CD4, integrated_cluster.IDs == "CD4+ T" )
DefaultAssay(NonnaiveCD4) <- "RNA"
NonnaiveCD4.postBoost <- subset(NonnaiveCD4, timepoint == "1moPostBoost" & sample == "unstim")

Idents(NonnaiveCD4.postBoost)<- "PriorCOVID"
table(NonnaiveCD4.postBoost$PriorCOVID)

NonnaiveCD4.postBoost.unstim <- FindMarkers(NonnaiveCD4.postBoost, ident.1 = "yes", ident.2 = "no", test.use="DESeq2") 
NonnaiveCD4.postBoost.unstim$diffexpressed <- "NO"
NonnaiveCD4.postBoost.unstim$diffexpressed[NonnaiveCD4.postBoost.unstim$avg_log2FC > 0.05 & 
                                             NonnaiveCD4.postBoost.unstim$p_val_adj < 0.05] <- "UP"
NonnaiveCD4.postBoost.unstim$diffexpressed[NonnaiveCD4.postBoost.unstim$avg_log2FC < -0.05 & 
                                             NonnaiveCD4.postBoost.unstim$p_val_adj < 0.05] <- "DOWN"

#write.csv(NonnaiveCD4.postBoost.unstim,"./Gene_sets/Unstim_NonnaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = TRUE)
NonnaiveCD4.postBoost.unstim <- read.csv("./Gene_sets/Unstim_NonnaiveCD4_infection_UP_vax_DOWN_PostBoost.DESeq.csv", row.names = 1)

Experienced.unstimCD4.DEseq <- subset(NonnaiveCD4.postBoost.unstim, subset = diffexpressed == "UP" )
Naive.unstimCD4.DEseq <- subset(NonnaiveCD4.postBoost.unstim, subset = diffexpressed == "DOWN" )

top <- rbind(Experienced.unstimCD4.DEseq,Naive.unstimCD4.DEseq)%>% row.names()
plot2<- ggplot(data=NonnaiveCD4.postBoost.unstim, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 1)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,2)) +scale_y_continuous(limits = c(-0.1,100))+ 
  scale_color_manual(values = c(pre_vax,"lightgrey",pre_infect))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 6, max.overlaps = 15)+labs(title="Post Boost")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z

top <- rbind(Experienced.unstimCD4.DEseq,Naive.unstimCD4.DEseq)
z<- ggplot(top, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(pre_vax,pre_infect))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=after_stat(count)), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,670)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/Unstim_NonnaiveCD4_PostBoost_PriorCOVID_DEGsummary_BarGraph.jpeg", width = 6, height =6)

infection_list <- list(ISG,  prolif, 
                       rownames(Infection_imprint), rownames(Vaccine_imprint),
                       rownames(Experienced.unstimCD4.DEseq), rownames(Naive.unstimCD4.DEseq))
names(infection_list) = c("ISG","Prolif",
                          "Infection_imprint", "Vaccine_imprint", 
                          "unstim_infection_gex_up","unstim_vaccine_gex_up")

# Pairwise venn diagram
ggvenn(infection_list, c("Infection_imprint","unstim_infection_gex_up"),fill_color = c(infect, "lightblue"),auto_scale = T, 
       set_name_size = 0)          
ggvenn(infection_list, c("Vaccine_imprint","unstim_vaccine_gex_up"),fill_color = c( vax, "lightblue"), auto_scale = T, 
       set_name_size = 0)




#---------HEATMAP OF differentially expressed genes between experienced and naive in CD4 AIM reactive cells post boost-------------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint == "1moPostBoost" )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ PriorCOVID) 
#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","1moPostBoost_yes"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_vax, post_infect )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~PriorCOVID, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
z<- plotPCA(vsd, intgroup = "subgroup")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + scale_color_manual(values = c(post_vax, post_infect )) #"#ffa527","#56e39d", "#716be4" )) 
z
ggsave(z, filename = "~/S3data/Images/Images/AIMreactiveCD4_PriorCOVID_PCA.jpeg",width = 4, height= 4)

pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
top <- rbind(Experienced.AIMreactive.CD4.DEseq,Naive.AIMreactive.CD4.DEseq)%>% row.names()
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "PriorCOVID")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  PriorCOVID = c(yes = post_infect, no = post_vax ))
div_colors <- hcl_palettes(type = "diverging",palette = "Blue-Red 3")
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA,
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        #show_rownames = T,
                        #labels_row = genes,
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        annotation_legend = F,
                        annotation_names_col = FALSE, 
                        color = colorRampPalette(c("#002F70", "white", "#881415"))(50))#

z

genes = c("GZMB","CCL5","CCL4","HLA-B","GZMH","NKG7","IFITM1","CCL4","IFITM3","CCL3","GZMB","IFI6","LTB","CST7",
          "SERPINB9","CCR7","NFKB2","DDIT4","NFKBIA","NINJ1","GRAMD1A","NFKBIZ","ICOS","IL2","NFKBID","REL","STAT3")

z<- add.flag(z,
         kept.labels = genes,
         repel.degree = .75)

ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_DEGs_PostBoost.pdf", width = 4, height = 5)

#---------Violin plots for postbreakthru--------------
AIM.reactive.CD4.breakthru <- subset(AIM.reactive.CD4, timepoint != "7moPostVax")
AIM.reactive.CD4.breakthru@meta.data <- AIM.reactive.CD4.breakthru@meta.data %>%
  unite("COVIDtimepoint", PriorCOVID,timepoint, remove = F)
DefaultAssay(AIM.reactive.CD4.breakthru) <- "RNA"
z <- VlnPlot(AIM.reactive.CD4.breakthru, features = c("GZMB","REL"),ncol = 1,
             group.by  = "COVIDtimepoint", pt.size = 0.0) & 
  scale_fill_manual(values = c(post_vax,breakthru,post_infect))# & labs(x = "", y = "") & 
  #theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_IFI6_breakthru.pdf", width = 3, height =2 )


AIM.reactive.CD4.breakthru <-subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" & PriorCOVID == "yes" | 
                                      timepoint == "PostBreakThru"  )
DefaultAssay(AIM.reactive.CD4.breakthru) <- "RNA"
AIM.reactive.CD4.breakthru$PriorCOVID <- factor(AIM.reactive.CD4.breakthru$PriorCOVID, levels = c("yes","no"))
z <- VlnPlot(AIM.reactive.CD4.breakthru, features = c("GZMB"),ncol = 1,
             group.by  = "PriorCOVID", pt.size = 0.0, y.max = 6) & 
  scale_fill_manual(values = c(post_infect,breakthru)) & NoLegend() & labs(x = "", y = "") & 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
z
ggsave(z, filename = "./Images/Images/VlnPlot_DEGs_GZMB_breakthru.pdf", width = 3, height =2 )

#---------find differentially expressed genes between naive post boost and breakthru in CD4 AIM reactive cells-------------
Idents(AIM.reactive.CD4)<- "timepoint"
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.breakthru <- subset(AIM.reactive.CD4, PriorCOVID == "no" & timepoint != "7moPostVax")
table(AIM.reactive.CD4.breakthru$timepoint,AIM.reactive.CD4.breakthru$PriorCOVID )
PriorCOVID.CD4 <- FindMarkers(AIM.reactive.CD4.breakthru, ident.1 = "PostBreakThru", ident.2 = "1moPostBoost", test.use="DESeq2") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"

PostBreakThru.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
PreBreakThru.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq <- PriorCOVID.CD4

top <- rbind(PostBreakThru.AIMreactive.CD4.DEseq,PreBreakThru.AIMreactive.CD4.DEseq)%>% row.names()
top <- c("IFITM1","IFI6","RPS29")
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 2)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,1.5)) +scale_y_continuous(limits = c(-0.1,125))+ 
  scale_color_manual(values = c(post_vax,"lightgrey",breakthru))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0.1, ynudge = 0.1, size = 6)+labs(title="Post Breakthru")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_VaxPrimed_Volcano.pdf", width = 6, height =6)
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_VaxPrimed_Volcano.jpeg", width = 6, height =6)
#write.csv(PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq,"./Gene_sets/PostBreakThru_UP_PreBreakThru_DOWN_CD4_VaccinePrimed.DESeq.csv", row.names = TRUE)
PriorCOVID.CD4 <- read.csv("./Gene_sets/PostBreakThru_UP_PreBreakThru_DOWN_CD4_VaccinePrimed.DESeq.csv", row.names = 1 )

top.1 <- rbind(PostBreakThru.AIMreactive.CD4.DEseq,PreBreakThru.AIMreactive.CD4.DEseq)
z<- ggplot(top.1, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(post_vax,breakthru))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=..count..), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,Subject30)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_VaxPrimed_DEGsummary_BarGraph.jpeg", width = 6, height =6)
#---------HEATMAP OF differentially expressed genes CD4 AIM reactive cells post breakthru-------------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.sce <- AIM.reactive.CD4.breakthru %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","PostBreakThru_no"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_vax, breakthru  )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~timepoint, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
z<- plotPCA(vsd, intgroup = "subgroup")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + 
  scale_color_manual(values = c(post_vax, breakthru))  
z

top <- rbind(PostBreakThru.AIMreactive.CD4.DEseq,PreBreakThru.AIMreactive.CD4.DEseq)%>% row.names()
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "subgroup")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  subgroup = c(`1moPostBoost_no` = post_vax, `PostBreakThru_no` = breakthru))


z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,cellheight = 0.5,
                        border_color = NA,cutree_rows = 1,
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        #show_rownames = T,
                        #labels_row = genes,
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        annotation_legend = F,
                        annotation_names_col = FALSE,
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))#

z

genes = c("CCR7","IFITM1","IFI6","IFITM3","ISG15","ISG20","NFKBIA")

z<- add.flag(z,
             kept.labels = genes,
             repel.degree = 0)

ggsave(z, filename = "./Images/Images/Heatmap_DEGs_VAX_vs_PostBreakthrut.pdf", width = 3.5, height = 2)                                       
ggsave(z, filename = "./Images/Images/Heatmap_DEGs_VAX_vs_PostBreakthrut.jpeg", width = 3.5, height = 2)                                       


#---------find differentially expressed genes between infection post boost and breakthru in CD4 AIM reactive cells-------------
Idents(AIM.reactive.CD4)<- "timepoint"
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.breakthru <-subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" & PriorCOVID == "yes" | 
                                                              timepoint == "PostBreakThru"  )
table(AIM.reactive.CD4.breakthru$timepoint,AIM.reactive.CD4.breakthru$PriorCOVID )
PriorCOVID.CD4 <- FindMarkers(AIM.reactive.CD4.breakthru, ident.1 = "PostBreakThru", ident.2 = "1moPostBoost", test.use="DESeq2") 
PriorCOVID.CD4$diffexpressed <- "NO"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC > 0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "UP"
PriorCOVID.CD4$diffexpressed[PriorCOVID.CD4$avg_log2FC < -0.05 & 
                               PriorCOVID.CD4$p_val_adj < 0.05] <- "DOWN"

PostBreakThru.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "UP" )
PreBreakThru.AIMreactive.CD4.DEseq <- subset(PriorCOVID.CD4, subset = diffexpressed == "DOWN" )
PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq <- PriorCOVID.CD4

top <- rbind(PostBreakThru.AIMreactive.CD4.DEseq,PreBreakThru.AIMreactive.CD4.DEseq)%>% row.names()
top.left <- c("SNHG5","RPS26", "HLA-B","HLA-A","CCL5","GZMH","GZMB","CST7","CCL4","CCL3")
top.right <- c("DDIT4","NINJ1","TRAF4", "CCR7","INSIG1","IFITM2")
plot2<- ggplot(data=PriorCOVID.CD4, aes(x=avg_log2FC, y=-log10(p_val_adj), col =diffexpressed))+
  #geom_hline(yintercept = 1.30103, linetype = "longdash", color = "lightgrey")+
  #geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point(size = 2)+ theme_classic()+ scale_x_continuous(limits = c(-1.5,1.5)) + #scale_y_continuous(limits = c(-0.1,125))+ 
  scale_color_manual(values = c(post_infect,"lightgrey",breakthru))
z <-LabelPoints(plot = plot2, points = top.left, repel = TRUE,xnudge =-0.1, ynudge=0,size = 6)+labs(title="Post Breakthru")+
  NoLegend()+labs(y = "-log10 adjusted p value", x = "log2 (Fold Change)")+
  theme() + theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
                  axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
                  legend.title =element_text(size = 15), axis.title.x = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) 
z <- LabelPoints(plot = z, 
              points = top.right, repel = TRUE,xnudge =0.1, ynudge=0,size = 6)
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_InfectPrimed_Volcano.pdf", width = 6, height =6)
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_InfectPrimed_Volcano.jpeg", width = 6, height =6)
#write.csv(PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq,"./Gene_sets/PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq.csv", row.names = TRUE)
PriorCOVID.CD4 <- read.csv("./Gene_sets/PostBreakThru_UP_PreBreakThru_DOWN_CD4.DESeq.csv", row.names = 1)
View(PriorCOVID.CD4)

top.2 <- rbind(PostBreakThru.AIMreactive.CD4.DEseq,PreBreakThru.AIMreactive.CD4.DEseq)
z<- ggplot(top.2, aes(x=diffexpressed))+geom_bar(stat="count", aes(fill = diffexpressed)) + scale_fill_manual(values = c(post_infect,breakthru))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=20)) +labs(y = "Number of DEGs", x ="") + 
  geom_text(stat='count',aes(label=..count..), vjust = -0.5, size = 15) + scale_y_continuous(limits=c(0,Subject30)) +NoLegend()
z
ggsave(z, filename = "./Images/Images/AIMreactive_CD4_PostBreakThru_InfectPrimed_DEGsummary_BarGraph.jpeg", width = 6, height =6)

x<- c(44,42,202,614) %>%  as.data.frame()
x$sample <- c("Vaccine-primed","Vaccine-primed","Infection-primed","Infection-primed")
x$diffex <- c("vax","breakthru","infect","breakthru") 
x$diffex <- factor(x$diffex, levels=c("breakthru","vax","infect"))
colnames(x) <- c("DEGs","sample","diffex")

z<- ggplot(x, aes(x=sample, y = DEGs))+geom_bar(stat="identity", aes(fill = sample)) + 
  scale_fill_manual(values = c(post_infect,post_vax))+
  theme_classic()+theme(axis.text=element_text(colour="black",size=20),
                        axis.title.y =element_text(colour="black",size=15)) +labs(y = "Number of DEGs", x ="") +NoLegend()+
  scale_y_continuous(limits=c(0,1000))
z
ggsave(z, filename = "./Images/Images/DEGs_BarGraph_Figure5.pdf", width = 4, height = 6)                                            

#---------HEATMAP OF differentially expressed genes CD4 AIM reactive cells post breakthru-------------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.sce <- AIM.reactive.CD4.breakthru %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ PriorCOVID) 
#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_yes","PostBreakThru_no"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_infect, breakthru  )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~PriorCOVID, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
z<- plotPCA(vsd, intgroup = "subgroup")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + scale_color_manual(values = c(post_infect, breakthru))  
z

top <- rbind(top.2) %>% row.names()
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "subgroup")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  subgroup = c(`1moPostBoost_no` = post_vax, `PostBreakThru_no` = breakthru, `1moPostBoost_yes` = post_infect ))

z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA, cellheight = 0.5,
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        #show_rownames = T,
                        #labels_row = genes,
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        annotation_legend = F,
                        annotation_names_col = FALSE,
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))#

z

genes = c("GZMB","CCL5","CCL4","REL","HLA-B","GZMH","NKG7","IFITM1","CCL4","IFITM3","CCL3","GZMB","IFI6","LTB","CST7",
          "SERPINB9","CCR7","NFKB2","DDIT4","NFKBIA","NINJ1","GRAMD1A","NFKBIZ","ICOS")

z<- add.flag(z,
             kept.labels = genes,
             repel.degree = .75)

ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_DEGs_PostBreakthrut.pdf", width = 3.5, height = 6)                                       
ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_DEGs_PostBreakthrut.jpeg", width = 3.5, height = 6)                                       

#---------Metascape analysis of Infection-primed CD4------------------
metasacpe_result <- read_excel("~/S3data/Metascape_analyses/METAscape_Experienced_enrichedTerms_PostBoost/metascape_result.xlsx", 
                               sheet = "Enrichment")
subset <- metasacpe_result %>%
  filter(LogP <= -1.3, str_ends(GroupID, "Summary"))
View(subset)
subset <- subset %>%
  filter(Description %in% c("Cytokine Signaling in Immune system",
                            "Interferon Signaling",
                            "Natural killer cell mediated cytotoxicity",
                            "Type II interferon signaling",
                            "inflammatory response", "Network map of SARS-CoV-2 signaling pathway", "response to interferon-beta"))

z<- ggplot(subset, aes(x=-LogP, y=reorder(Description, -LogP), fill = -LogP))+geom_col()+ 
  scale_fill_gradient(low = "lightgrey",high = "purple4",limits=c(0,30))+
  theme(legend.key = element_rect(fill=NA)) + 
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1),
        axis.text=element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15))+labs(y="",x="")
z
ggsave(z, filename = "./Images/Images/EnrichedTerms_Experienced_PostBoost.pdf", width = 7, height =3)


metasacpe_result <- read_excel("~/S3data/Metascape_analyses/METAscape_Naive_enrichedTerms_PostBoost/metascape_result.xlsx", 
                               sheet = "Enrichment")
#View(metasacpe_result)
#subset <- metasacpe_result %>%
#  filter(LogP <= -6.75, str_ends(GroupID, "Summary"))
subset <- metasacpe_result %>%
  filter(Description %in% c("Cell Cycle, Mitotic",
                            "Cytokine Signaling in Immune system",
                            "Mitotic G2-G2/M phases",
                            "Signaling by Interleukins",
                            "G2/M Checkpoints", "regulation of I-kappaB kinase/NF-kappaB signaling", "T cell activation"))


z<- ggplot(subset, aes(x=-LogP, y=reorder(Description, -LogP), fill = -LogP))+geom_col()+ 
  scale_fill_gradient(low = "lightgrey",high = "orange3",limits=c(0,28))+
  theme(legend.key = element_rect(fill=NA)) + 
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1),
        axis.text=element_text(colour="black",size=15), legend.text=element_text(colour="black",size=15))+labs(y="",x="")
z
ggsave(z, filename = "./Images/Images/EnrichedTerms_Vaccine_PostBoost.pdf", width = 7.25, height =3)

#---------pseudo bulk of AIM reactive CD4 by subject at PreBoost---------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint == "7moPostVax" )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ PriorCOVID) 
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("7moPostVax_no","7moPostVax_yes")) #"1moPostBoost_no","PostBreakThru_no","1moPostBoost_yes"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c("#ffc471","#b5b2f1"))#"#ffa527","#56e39d", "#716be4" )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
 plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values = c("#ffc471","#b5b2f1"))#"#ffa527","#56e39d", "#716be4" )) 
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 3 # filter for genes with at least 10 counts in 3/14 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("PriorCOVID","yes","no"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")]
write.table(geneRank, file = "~/S3data/Gene_sets/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F,sep = "\t")

#volcano
View(diffExpr)
diffExpr$diffexpressed <- "NO"
diffExpr$diffexpressed[diffExpr$log2FoldChange > 0.1 & 
                         diffExpr$padj <= 0.05] <- "UP"
diffExpr$diffexpressed[diffExpr$log2FoldChange < -0.1 & 
                         diffExpr$padj <= 0.05] <- "DOWN"
UP <- subset(diffExpr, subset = diffexpressed == "UP" )
DOWN <- subset(diffExpr, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
plot2<- ggplot(data=diffExpr, aes(x=log2FoldChange, y=-log10(padj), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+ geom_vline(xintercept = c(-0.1, 0.1), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+  scale_color_manual(values = c("lightgrey","#b5b2f1"))
plot2
#ggsave(z, filename = "~/S3data/Images/Images/PreBoost_PseudoBulk_volcano.pdf",width = 8, height = 8)
#ggsave(z, filename = "~/S3data/Images/Images/PreBoost_PseudoBulk_volcano.jpeg",width = 8, height = 8)

View(diffExpr)
diffExpr<- rownames_to_column(diffExpr, var = "genes")
top.yes <-  diffExpr %>% arrange((stat)) %>% head(10)
top.no <-  diffExpr %>% arrange(-(stat)) %>% head(10)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix
metaData <- metaData[,c("subject_timepoint", "PriorCOVID")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  PriorCOVID = c(no = "#ffc471",yes = "#b5b2f1" ))
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA,
                        annotation_col = metaData, annotation_colors	= cols,  # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        fontsize = 10, # Shrink the pathway labels
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

z
ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_rankedDEseq_PreBoost.pdf", width = 4, height = 5)
ggsave(z, filename = "./Images/Images//Heatmap_PriorCOVID_rankedDEseq_PreBoost.jpeg", width = 4, height = 5)

#---------pseudo bulk of AIM reactive CD4 by subject at PostBoost---------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint == "1moPostBoost" )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ PriorCOVID) 
#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","1moPostBoost_yes"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c("#ffa527", "#716be4" )) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~PriorCOVID, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + 
  scale_color_manual(values = c("#ffa527", "#716be4" )) #"#ffa527","#56e39d", "#716be4" )) 
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 3 # filter for genes with at least 10 counts in 3/14 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("PriorCOVID","yes","no"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")]
View(geneRank)
write.table(geneRank, file = "~/S3data/Gene_sets/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F, sep = "\t")

#volcano
View(diffExpr)
diffExpr$diffexpressed <- "NO"
diffExpr$diffexpressed[diffExpr$log2FoldChange > 0.1 & 
                         diffExpr$padj <= 0.05] <- "UP"
diffExpr$diffexpressed[diffExpr$log2FoldChange < -0.1 & 
                         diffExpr$padj <= 0.05] <- "DOWN"
UP <- subset(diffExpr, subset = diffexpressed == "UP" )
DOWN <- subset(diffExpr, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
plot2<- ggplot(data=diffExpr, aes(x=log2FoldChange, y=-log10(padj), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+ geom_vline(xintercept = c(-0.1, 0.1), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+  scale_color_manual(values = c("#ffa527","lightgrey","#716be4"))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 4)+labs(title="AIM Reactive CD4+ T cells (pseudobulk, post boost)")+
  NoLegend()+labs(y = "-log10 p value", x = "Naive -------------- log2 (Fold Change) -------------- Experienced")+
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1))
z
ggsave(z, filename = "~/S3data/Images/Images/PostBoost_PseudoBulk_volcano.pdf",width = 8, height = 8)
ggsave(z, filename = "~/S3data/Images/Images/PostBoost_PseudoBulk_volcano.jpeg",width = 8, height = 8)

View(diffExpr)
diffExpr<- rownames_to_column(diffExpr, var = "genes")
top.yes <-  diffExpr %>% arrange((stat)) %>% head(10)
top.no <-  diffExpr %>% arrange(-(stat)) %>% head(10)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "PriorCOVID")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  PriorCOVID = c(yes = "#716be4", no = "#ffa527" ))

z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        cutree_row = 2,
                        border_color = NA,
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        #show_rownames = T,
                        #labels_row = genes,
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        annotation_legend = F,
                        annotation_names_col = FALSE,
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))#
z

ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_rankedDEseq_PostBoost.pdf", width = 4, height = 5)
ggsave(z, filename = "./Images/Images//Heatmap_PriorCOVID_rankedDEseq_PostBoost.jpeg", width = 4, height = 5)
#---------pseudo bulk of AIM reactive CD4 by subject at PostBreakthru---------
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" & PriorCOVID == "no"  )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ subject + timepoint) 
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","PostBreakThru_no"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_vax,breakthru)) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values = c(post_vax,breakthru))
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 3 # filter for genes with at least 10 counts in 2/12 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("timepoint","PostBreakThru","1moPostBoost"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")] %>% as.data.frame()
View(geneRank)
write.table(geneRank, file = "~/S3data/Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
             col.names = F, quote =F, sep = "\t")

#volcano
View(diffExpr)
diffExpr$diffexpressed <- "NO"
diffExpr$diffexpressed[diffExpr$log2FoldChange > 0.1 & 
                         diffExpr$padj <= 0.05] <- "UP"
diffExpr$diffexpressed[diffExpr$log2FoldChange < -0.1 & 
                         diffExpr$padj <= 0.05] <- "DOWN"
UP <- subset(diffExpr, subset = diffexpressed == "UP" )
DOWN <- subset(diffExpr, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
plot2<- ggplot(data=diffExpr, aes(x=log2FoldChange, y=-log10(padj), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+ geom_vline(xintercept = c(-0.1, 0.1), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+  scale_color_manual(values = c("lightgrey"))
plot2
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.pdf",width = 8, height = 8)
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.jpeg",width = 8, height = 8)

View(diffExpr)
diffExpr<- rownames_to_column(diffExpr, var = "genes")
top.yes <-  diffExpr %>% arrange((stat)) %>% head(10)
top.no <-  diffExpr %>% arrange(-(stat)) %>% head(10)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix
metaData <- metaData[,c("subject_timepoint", "timepoint")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  timepoint = c(PostBreakThru = "#56e39d" ,`1moPostBoost` = "#ffc471"))
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA,
                        annotation_col = metaData, annotation_colors	= cols,  # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        fontsize = 10, # Shrink the pathway labels
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

z

#---------pseudo bulk of AIM reactive CD4 by subject at PostBreakthru ONLY with paired samples---------
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" & PriorCOVID == "no"  )
table(AIM.reactive.CD4.sce$subject, AIM.reactive.CD4.sce$timepoint)
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4.sce, subject == "Subject1" |
                                 subject == "Subject2" | subject == "Subject4" | 
                                 subject == "Subject5"| subject == "Subject6")
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ subject + timepoint) 
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","PostBreakThru_no"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c("#ffc471","#56e39d")) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values = c("#ffc471","#56e39d"))
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 2 # filter for genes with at least 10 counts in 2/10 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("timepoint","PostBreakThru","1moPostBoost"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")] %>% as.data.frame()
View(geneRank)
write.table(geneRank, file = "~/S3data/Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq.JUSTBREAKTHRUS.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F, sep = "\t")

#volcano
View(diffExpr)
diffExpr$diffexpressed <- "NO"
diffExpr$diffexpressed[diffExpr$log2FoldChange > 0.1 & 
                         diffExpr$padj <= 0.05] <- "UP"
diffExpr$diffexpressed[diffExpr$log2FoldChange < -0.1 & 
                         diffExpr$padj <= 0.05] <- "DOWN"
UP <- subset(diffExpr, subset = diffexpressed == "UP" )
DOWN <- subset(diffExpr, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
plot2<- ggplot(data=diffExpr, aes(x=log2FoldChange, y=-log10(padj), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+ geom_vline(xintercept = c(-0.1, 0.1), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+  scale_color_manual(values = c("lightgrey"))
plot2
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.pdf",width = 8, height = 8)
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.jpeg",width = 8, height = 8)

View(diffExpr)
diffExpr<- rownames_to_column(diffExpr, var = "genes")
top.yes <-  diffExpr %>% arrange((stat)) %>% head(10)
top.no <-  diffExpr %>% arrange(-(stat)) %>% head(10)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix
metaData <- metaData[,c("subject_timepoint", "timepoint")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  timepoint = c(PostBreakThru = "#56e39d" ,`1moPostBoost` = "#ffc471"))
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA,
                        annotation_col = metaData, annotation_colors	= cols,  # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        fontsize = 10, # Shrink the pathway labels
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

z


#---------pseudo bulk of AIM reactive CD4 by subject at PostBreakthru versus infection---------
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" & PriorCOVID == "yes" | 
                                 timepoint == "PostBreakThru"  )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_yes","PostBreakThru_no"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_infect,breakthru)) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values =  c(post_infect,breakthru)) 
#filter 
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) > 3 # filter for genes with at least 10 counts in 2/12 samples 
fullDataset <- dds[keep,]
DESdata_mixed <- DESeq(fullDataset, parallel=TRUE)

diffExpr <- results(DESdata_mixed, contrast=c("timepoint","PostBreakThru","1moPostBoost"), parallel = T) %>% as.data.frame()
write.csv(diffExpr,"./Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Experienced_PseudoBulk_DESeq.csv", row.names = TRUE)
geneRank <- rownames_to_column(diffExpr, "genes")
geneRank <- geneRank %>% arrange(-stat)
geneRank <- geneRank[,c("genes","stat")] %>% as.data.frame()
View(geneRank)
write.table(geneRank, file = "~/S3data/Gene_sets/BreakThru_UP_PostBoost_DOWN_CD4_Experienced_PseudoBulk_DESeq.rnk.txt",row.names = FALSE, 
            col.names = F, quote =F, sep = "\t")

#volcano
View(diffExpr)
diffExpr$diffexpressed <- "NO"
diffExpr$diffexpressed[diffExpr$log2FoldChange > 0.1 & 
                         diffExpr$padj <= 0.05] <- "UP"
diffExpr$diffexpressed[diffExpr$log2FoldChange < -0.1 & 
                         diffExpr$padj <= 0.05] <- "DOWN"
UP <- subset(diffExpr, subset = diffexpressed == "UP" )
DOWN <- subset(diffExpr, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
plot2<- ggplot(data=diffExpr, aes(x=log2FoldChange, y=-log10(padj), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+ geom_vline(xintercept = c(-0.1, 0.1), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+  scale_color_manual(values = c("lightgrey"))
plot2
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.pdf",width = 8, height = 8)
#ggsave(z, filename = "~/S3data/Images/Images/BreakThru_PseudoBulk_volcano.jpeg",width = 8, height = 8)

View(diffExpr)
diffExpr<- rownames_to_column(diffExpr, var = "genes")
top.yes <-  diffExpr %>% arrange((stat)) %>% head(10)
top.no <-  diffExpr %>% arrange(-(stat)) %>% head(10)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix
metaData <- metaData[,c("subject_timepoint", "timepoint")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  timepoint = c(PostBreakThru = "#56e39d" ,`1moPostBoost` = "#716be4"))
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = F,
                        border_color = NA,
                        annotation_col = metaData, annotation_colors	= cols,  # Add metadata labels!
                        show_colnames = F, # Don't show sample labels
                        fontsize = 10, # Shrink the pathway labels
                        color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

z
ggsave(z, filename = "./Images/Images/Heatmap_Breakthru_infectionPrimed_rankedDEseq.pdf", width = 4, height = 5)

#---------peudo bulk of AIM reactive CD4 by subject at PostBoost and PostBreakththru for sample PCA------------
AIM.reactive.CD4.sce <- subset(AIM.reactive.CD4, sample == "stim" & timepoint != "7moPostVax" )
AIM.reactive.CD4.sce <- AIM.reactive.CD4.sce %>% as.SingleCellExperiment() #using RNA counts !!!

## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")],use.assay.type = "counts")
#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes", endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN", endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",endsWith(subject_timepoint, "Subject2_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",endsWith(subject_timepoint, "Subject6_PostBreakThru")  ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH", endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject_timepoint, "Subject5_PostBreakThru")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ subject + timepoint) 
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=TRUE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("1moPostBoost_no","PostBreakThru_no", "1moPostBoost_yes"))
plotPCA(vsd, intgroup = "runDate")+ geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c(post_vax,breakthru, post_infect)) 
#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "runDate")+geom_point( size = 6)+ theme_classic()  
plotPCA(vsd, intgroup = "subgroup")+ geom_point( size = 2)+ theme_classic() + scale_color_manual(values = c(post_vax,breakthru, post_infect)) 
z<- plotPCA(vsd, intgroup = "subgroup")+ geom_point(size = 3, color = "black", shape = 21)+ theme_classic() + 
  scale_color_manual(values = c(post_vax,breakthru, post_infect))  #"#ffa527","#56e39d", "#716be4" )) 
 # scale_x_continuous(limits = c(-25,25)) + scale_y_continuous(limits = c(-10,10))
z

ggsave(z, filename = "./Images/Images/AIMreactiveCD4_BreakThru_PCA.jpeg", width = 5, height= 5)



#-------------------------------GSEA results-------------------------------
Experienced_UP_Naive_DOWN_CD4_PreBoost_pos <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/gsea_report_for_na_pos_1671501202408.tsv")
Experienced_UP_Naive_DOWN_CD4_PreBoost_neg <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/gsea_report_for_na_neg_1671501202408.tsv")
GSEA_PreBoost <- rbind(Experienced_UP_Naive_DOWN_CD4_PreBoost_pos,Experienced_UP_Naive_DOWN_CD4_PreBoost_neg)
GSEA_PreBoost$NAME<-gsub("HALLMARK_","",as.character(GSEA_PreBoost$NAME))
GSEA_PreBoost$NAME<- str_to_title(GSEA_PreBoost$NAME) 
View(GSEA_PreBoost)

Experienced_UP_Naive_DOWN_CD4_PostBoost_pos <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/gsea_report_for_na_pos_1671501213901.tsv")
Experienced_UP_Naive_DOWN_CD4_PostBoost_neg <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/gsea_report_for_na_neg_1671501213901.tsv")
GSEA_PostBoost <- rbind(Experienced_UP_Naive_DOWN_CD4_PostBoost_neg,Experienced_UP_Naive_DOWN_CD4_PostBoost_pos)
GSEA_PostBoost$NAME<-gsub("HALLMARK_","",as.character(GSEA_PostBoost$NAME))
GSEA_PostBoost$NAME<- str_to_title(GSEA_PostBoost$NAME) 

BreakTHRU_UP_Postboost_DOWN_CD4_Naive_pos <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq/gsea_report_for_na_pos_1671501220733.tsv")
BreakTHRU_UP_Postboost_DOWN_CD4_Naive_neg <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Naive_PseudoBulk_DESeq/gsea_report_for_na_neg_1671501220733.tsv")
GSEA_PostBreakthru <- rbind(BreakTHRU_UP_Postboost_DOWN_CD4_Naive_pos,BreakTHRU_UP_Postboost_DOWN_CD4_Naive_neg)
GSEA_PostBreakthru$NAME<- gsub("HALLMARK_","",as.character(GSEA_PostBreakthru$NAME)) 
GSEA_PostBreakthru$NAME<- str_to_title(GSEA_PostBreakthru$NAME) 


BreakTHRU_UP_Postboost_DOWN_CD4_Infect_pos <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Infect_PseudoBulk_DESeq/gsea_report_for_na_pos_1671501228634.tsv")
BreakTHRU_UP_Postboost_DOWN_CD4_Infect_neg <- 
  read_tsv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Infect_PseudoBulk_DESeq/gsea_report_for_na_neg_1671501228634.tsv")
GSEA_PostBreakthru_infect <- rbind(BreakTHRU_UP_Postboost_DOWN_CD4_Infect_pos,BreakTHRU_UP_Postboost_DOWN_CD4_Infect_neg)
GSEA_PostBreakthru_infect$NAME<- gsub("HALLMARK_","",as.character(GSEA_PostBreakthru_infect$NAME)) 
GSEA_PostBreakthru_infect$NAME<- str_to_title(GSEA_PostBreakthru_infect$NAME) 

Experienced_UP_Naive_DOWN_CD4_PreBoost.sub <- subset(GSEA_PreBoost, 
                                                     `FDR q-val` <= 0.05 )
z<-ggplot(Experienced_UP_Naive_DOWN_CD4_PreBoost.sub,aes(x=NES,y= reorder(NAME,-`NES`))) + 
  geom_point(aes(size=`FDR q-val`)) + scale_size(range = c(1, 5), name="FDR q-val", trans = 'reverse') + 
  geom_segment(aes(xend=0, yend = NAME)) +
  theme(legend.key = element_rect(fill=NA)) + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1), 
        axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Pre Booster", x ="", 
       y = "")
z
ggsave(z, filename = "./Images/Images/GSEA_hallmark_Experienced-Naive_PreBooster.pdf",width = 7.25, height = 6)

Experienced_UP_Naive_DOWN_CD4_PostBoost.sub <- subset(GSEA_PostBoost, 
                                                      `FDR q-val` <= 0.05)
Experienced_UP_Naive_DOWN_CD4_PostBoost.sub <- subset(Experienced_UP_Naive_DOWN_CD4_PostBoost.sub, NAME != "E2f_targets" & NAME != "Estrogen_response_early"&
                                                        NAME !="Estrogen_response_late"&NAME != "Heme_metabolism"&NAME != "Allograft_rejection" )
z<-ggplot(Experienced_UP_Naive_DOWN_CD4_PostBoost.sub,aes(x=NES,y= reorder(NAME,-`NES`))) + 
  geom_point(aes(size=`FDR q-val`)) + scale_size(range = c(1, 5), name="FDR q-val", trans = 'reverse') + 
  geom_segment(aes(xend=0, yend = NAME)) +scale_x_continuous(limits = c(-2.5,3))+ 
  theme(legend.key = element_rect(fill=NA)) + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1), 
        axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Booster", x ="", 
       y = "")
z
ggsave(z, filename = "./Images/Images/GSEA_hallmark_Experienced-Naive_PostBooster.pdf", width = 7.25, height = 5.5)

GSEA_PostBreakthru.sub <- subset(GSEA_PostBreakthru, 
                                 `FDR q-val` <= 0.05)
GSEA_PostBreakthru.sub <- GSEA_PostBreakthru.sub %>%
  mutate(Color = case_when(
    NES > 0 ~ "breakthru",
    NES < 0 ~ "vax"
  ))
GSEA_PostBreakthru.sub$NAME
GSEA_PostBreakthru.sub <- subset(GSEA_PostBreakthru.sub, NAME != "Heme_metabolism"&NAME != "Allograft_rejection" )

z<-ggplot(GSEA_PostBreakthru.sub,aes(x=NES,y= reorder(NAME,-`NES`), color = Color, fill = Color)) + 
  geom_point(aes(size=`FDR q-val`)) + scale_size(range = c(1, 5), name="FDR q-val", trans = 'reverse') + 
  geom_segment(aes(xend=0, yend = NAME)) +scale_x_continuous(limits = c(-2.5,3)) +
  theme(legend.key = element_rect(fill=NA)) + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1),
        axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Breakthru", x ="", 
       y = "")+scale_color_manual(values = c(breakthru, post_vax))+scale_fill_manual(values = c(breakthru, post_vax))
z
ggsave(z, filename = "./Images/Images/GSEA_hallmark_BreakThru_Naive.pdf", width = 7, height = 5)

GSEA_PostBreakthru.sub <- subset(GSEA_PostBreakthru_infect, 
                                 `FDR q-val` <= 0.05)
GSEA_PostBreakthru.sub <- GSEA_PostBreakthru.sub %>%
  mutate(Color = case_when(
    NES > 0 ~ "breakthru",
    NES < 0 ~ "infect"
  ))
GSEA_PostBreakthru.sub$NAME
GSEA_PostBreakthru.sub <- subset(GSEA_PostBreakthru.sub, NAME != "E2f_targets" & NAME != "Estrogen_response_early"&
                                   NAME !="Estrogen_response_late"&NAME != "Heme_metabolism"&NAME != "Allograft_rejection" )

z<-ggplot(GSEA_PostBreakthru.sub,aes(x=NES,y= reorder(NAME,-`NES`), color = Color, fill = Color)) + 
  geom_point(aes(size=`FDR q-val`)) + scale_size(range = c(1, 5), name="FDR q-val", trans = 'reverse') + 
  geom_segment(aes(xend=0, yend = NAME)) + scale_x_continuous(limits = c(-4,3)) +
  theme(legend.key = element_rect(fill=NA)) + 
  theme(axis.text.x = element_text(size = 15), panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1), 
        axis.text.y = element_text(size = 20), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Breakthru", x ="", 
       y = "")+scale_color_manual(values = c(breakthru, post_infect))+scale_fill_manual(values = c(breakthru, post_infect))
z
ggsave(z, filename = "./Images/Images/GSEA_hallmark_BreakThru_Infect.pdf", width = 7.25, height = 5)

GSEA_PreBoost$timepoint <- "Pre"
#GSEA_PreBoost$NAME = paste0('Pre_', GSEA_PreBoost$NAME)
GSEA_PostBoost$timepoint <- "Post"
#GSEA_PostBoost$NAME = paste0('Post_', GSEA_PostBoost$NAME)

all <- rbind(GSEA_PostBoost,GSEA_PreBoost)
all$`FDR q-val`<- all$`FDR q-val` + 0.0001
View(all)
#all$`FDR q-val.x`<- all$`FDR q-val.x` + 0.001
#all$`FDR q-val.y`<- all$`FDR q-val.y` + 0.001
order <- GSEA_PostBoost %>% arrange(-(NES)) %>% as.data.frame()
all$timepoint <- factor(all$timepoint, levels = c( "Post", "Pre"))
all$NAME <- factor(all$NAME, levels = order$NAME)

z<-ggplot(all,aes(x=NES,y= NAME, color = timepoint, fill = timepoint)) + 
  scale_x_continuous(limits = c(-2.5,3)) +
  geom_bar(stat="identity", width =0.1, position = position_dodge(width = 0.5))+ 
  geom_point(aes(size=`FDR q-val`),position = position_dodge(width = 0.5)) + 
  scale_size(range = c(0,4), trans=c("log10","reverse"), name="FDR q-val") + 
  #geom_segment(aes(xend=0, yend = NAME, color = timepoint), position = position_dodge(width = -1))+
  theme(legend.key = element_rect(fill=NA),panel.background=element_rect(fill="white")) + 
  theme(axis.text.x = element_text(size = 15),
        panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed", size= 0.25, position_dodge(0.5)),
        panel.border=element_rect(fill = NA, colour = "black", size = 0.5), 
        axis.text.y = element_text(size = 10), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Breakthru", x ="", 
       y = "")+scale_color_manual(values = c("black", "grey53"))+scale_fill_manual(values = c("black","grey53"))
z
ggsave(z, filename = "./Images/Images/GSEA_hallmark_all.pdf", width = 6, height = 12)

View(all)

GSEA_PostBreakthru$test <- "vax"
GSEA_PostBreakthru_infect$test <- "infect"

all <- rbind(GSEA_PostBreakthru,GSEA_PostBreakthru_infect)
#all$`FDR q-val`<- all$`FDR q-val` + 0.0001
#all$`FDR q-val.x`<- all$`FDR q-val.x` + 0.001
#all$`FDR q-val.y`<- all$`FDR q-val.y` + 0.001
order <- GSEA_PostBreakthru_infect %>% arrange((NES)) %>% as.data.frame()
all$test <- factor(all$test, levels = c( "vax", "infect"))
all$NAME <- factor(all$NAME, levels = order$NAME)
all <- subset(all, NAME == "Hypoxia" | NAME == "Interferon_alpha_response"| NAME == "Interferon_gamma_response"| NAME =="Mitotic_spindle"|
                NAME  == "Glycolysis" | NAME == "Tgf_beta_signaling" |  NAME == "Tnfa_signaling_via_nfkb"| NAME == "Heme_metabolism"| 
                NAME == "Apoptosis"| NAME == "Allograft_rejection"| NAME == "Myc_targets_v1"| NAME == "Myc_targets_v2"| NAME == "Mtorc1_signaling"| 
                NAME == "Estrogen_response_late"| NAME == "Kras_signaling_dn"| NAME == "Estrogen_response_early")
all <- subset(all, `FDR q-val` <= 0.05)
z<-ggplot(all,aes(x=NES,y= NAME, color = test, fill = test)) + 
  scale_x_continuous(limits = c(-4,3)) +
  geom_bar(stat="identity", width =0.05, position = position_dodge(width = 0.5))+ 
  geom_point(aes(size=`FDR q-val`),position = position_dodge(width = 0.5)) + 
  scale_size(range = c(1,5), trans=c("reverse"), name="FDR q-val") + 
  #geom_segment(aes(xend=0, yend = NAME, color = timepoint), position = position_dodge(width = -1))+
  theme(legend.key = element_rect(fill=NA),panel.background=element_rect(fill="white")) + 
  theme(axis.text.x = element_text(size = 15),
        panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed", size= 0.25, position_dodge(0.5)),
        panel.border=element_rect(fill = NA, colour = "black", size = 0.5), 
        axis.text.y = element_text(size = 15), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Breakthru", x ="", 
       y = "")+scale_color_manual(values = c(post_vax,post_infect))+scale_fill_manual(values = c(post_vax,post_infect))+
  facet_grid(rows = "test")
z
View(all)
ggsave(z, filename = "./Images/Images/GSEA_hallmark_breakthru_summary.pdf", width = 6, height = 5)

all <- rbind(GSEA_PostBreakthru,GSEA_PostBreakthru_infect)
all$`FDR q-val`<- all$`FDR q-val` + 0.0001
#all$`FDR q-val.x`<- all$`FDR q-val.x` + 0.001
#all$`FDR q-val.y`<- all$`FDR q-val.y` + 0.001
order <- GSEA_PostBreakthru %>% arrange(-(NES)) %>% as.data.frame()
all$test <- factor(all$test, levels = c( "vax", "infect"))
all$NAME <- factor(all$NAME, levels = order$NAME)
z<-ggplot(all,aes(x=NES,y= NAME, color = test, fill = test)) + 
  scale_x_continuous(limits = c(-3.5,2.75)) +
  geom_bar(stat="identity", width =0.1, position = position_dodge(width = 0.5))+ 
  geom_point(aes(size=`FDR q-val`),position = position_dodge(width = 0.5)) + 
  scale_size(range = c(0.1,5), trans=c("log10","reverse"), name="FDR q-val") + 
  #geom_segment(aes(xend=0, yend = NAME, color = timepoint), position = position_dodge(width = -1))+
  theme(legend.key = element_rect(fill=NA),panel.background=element_rect(fill="white")) + 
  theme(axis.text.x = element_text(size = 15),
        panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed", size= 0.25, position_dodge(0.5)),
        panel.border=element_rect(fill = NA, colour = "black", size = 0.5), 
        axis.text.y = element_text(size = 10), legend.text =element_text(size = 15),
        legend.title =element_text(size = 15)) + geom_vline(xintercept = 0)+
  labs(title = "Post Breakthru", x ="", 
       y = "")+scale_color_manual(values = c(post_vax,post_infect))+scale_fill_manual(values = c(post_vax,post_infect))
ggsave(z, filename = "./Images/Images/GSEA_hallmark_breakthru_all.pdf", width = 6, height = 12)


#GSEA plots 
plotGSEAinfection <- function (Post,postPathway, pathwayName, title, leftLabel, rightLabel, cohortCompare, legendUpperRight)
{
  if (cohortCompare == F)     {     colorGrob1 <- "black"; colorGrob2 <- "black"     }
  if (cohortCompare == T)   { colorGrob1 <- "#ffa527"; colorGrob2 <- "#716be4"}
  annotationInfo <- paste0("NES: ", round(Post$NES[grep(pathwayName,Post$NAME)],2), "\n", "FDR: ", formatC(Post$`FDR q-val`[grep(pathwayName,Post$NAME)], format = "e", digits = 1))
  if (legendUpperRight == T)   {   grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.35, hjust=0, gp=gpar(col=post_infect, fontsize=20))) }
  if (legendUpperRight == F)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.25, hjust=0, gp=gpar(col="black", fontsize=24)))     }
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.97,hjust=0, gp=gpar(col=colorGrob1, fontsize=20)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.7,y=0.97,hjust=0, gp=gpar(col=colorGrob2, fontsize=20))) 
  return(
    ggplot(data=postPathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES) ) + geom_line(color=post_infect, size=1) + 
      geom_rug(sides="b", size=0.75, alpha=0.5, color=post_infect ) + theme_bw() +
      ggtitle(title) + ylab("Enrichment score") + xlab("Rank in gene list") + 
      theme(axis.text = element_text(size=20,hjust = 0.75, color="black"), axis.title = element_text(size=24,hjust = 0.5, color="black"), plot.title = element_text(size=32,hjust = 0.5))+
      annotation_custom(grob1) + geom_hline(yintercept = 0) + 
      annotation_custom(left_grob) + annotation_custom(right_grob) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
      geom_hline(yintercept = 0, size=0.5)
  )
}

IFNa_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Infect_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")
View(IFNa_post)
a<- plotGSEAinfection(GSEA_PostBreakthru_infect, IFNa_post, pathwayName = "Interferon_alpha_response", title = "Interferon Alpha Response", 
                        leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
a

IFNg_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/BreakThru_UP_PostBoost_DOWN_CD4_Infect_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")
b<- plotGSEAinfection(GSEA_PostBreakthru_infect,IFNg_post, pathwayName = "Interferon_gamma_response", title = "Interferon Gamma Response", 
                        leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
b

ggsave(a, filename = "./Images/Images/GSEA_IFNa_BreakThru.pdf", width = 7, height = 4.5)
ggsave(b, filename = "./Images/Images/GSEA_IFNg_BreakThru.pdf", width = 7, height = 4.5)



plotGSEAinfection <- function (Pre,Post, prePathway,postPathway, pathwayName, title, leftLabel, rightLabel, cohortCompare, legendUpperRight)
{
  if (cohortCompare == F)     {     colorGrob1 <- "black"; colorGrob2 <- "black"     }
  if (cohortCompare == T)   { colorGrob1 <- "#ffa527"; colorGrob2 <- "#716be4"}
  annotationInfo <- paste0("NES: ", round(Pre$NES[grep(pathwayName,Pre$NAME)],2), "\n", "FDR: ", formatC(Pre$`FDR q-val`[grep(pathwayName,Pre$NAME)], format = "e", digits = 1))
  annotationInfo2 <- paste0("NES: ", round(Post$NES[grep(pathwayName,Post$NAME)],2), "\n", "FDR: ", formatC(Post$`FDR q-val`[grep(pathwayName,Post$NAME)], format = "e", digits = 1))
  if (legendUpperRight == T)   {   grob1 = grobTree(textGrob(annotationInfo, x=0.68,  y=0.82, hjust=0, gp=gpar(col="grey53", fontsize=20)))
  grob2 = grobTree(textGrob(annotationInfo2, x=0.68,  y=0.52, hjust=0, gp=gpar(col="black", fontsize=20))) }
  if (legendUpperRight == F)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.25, hjust=0, gp=gpar(col="black", fontsize=24)))     }
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.97,hjust=0, gp=gpar(col=colorGrob1, fontsize=20)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.7,y=0.97,hjust=0, gp=gpar(col=colorGrob2, fontsize=20))) 
  return(
    ggplot(data=prePathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES) ) + geom_line(color="grey53", size=1) + 
      geom_rug(sides="t", size=0.75, alpha=0.5, color="grey53" ) + theme_bw() +
      geom_line(data=postPathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES),color="black", size=1)+ 
      geom_rug(data=postPathway, sides="b", size=0.75, alpha=0.5,color="black") +
      ggtitle(title) + ylab("Enrichment score") + xlab("Rank in gene list") + 
      theme(axis.text = element_text(size=19,hjust = 0.8, color="black"), axis.title = element_text(size=24,hjust = 0.5, color="black"), plot.title = element_text(size=32,hjust = 0.5))+
      annotation_custom(grob1) + annotation_custom(grob2)+ geom_hline(yintercept = 0) + 
      annotation_custom(left_grob) + annotation_custom(right_grob) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
      geom_hline(yintercept = 0, size=0.5)
  )
}

IFNa_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")
IFNa_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_ALPHA_RESPONSE.tsv", sep="\t")

a<- plotGSEAinfection(GSEA_PreBoost, GSEA_PostBoost, IFNa_pre,IFNa_post, pathwayName = "Interferon_alpha_response", title = "Interferon Alpha Response", 
                      leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
a

IFNg_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")
IFNg_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_INTERFERON_GAMMA_RESPONSE.tsv", sep="\t")

b<- plotGSEAinfection(GSEA_PreBoost, GSEA_PostBoost, IFNg_pre,IFNg_post, pathwayName = "Interferon_gamma_response", title = "Interferon Gamma Response", 
                      leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
b

Inflam_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_INFLAMMATORY_RESPONSE.tsv", sep="\t")
Inflam_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_INFLAMMATORY_RESPONSE.tsv", sep="\t")

e<- plotGSEAinfection(GSEA_PreBoost, GSEA_PostBoost, Inflam_pre,Inflam_post, pathwayName = "Inflammatory_response", title = "Inflammatory Response", 
                      leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
e

plotGSEAvaccine <- function (Pre,Post, prePathway,postPathway, pathwayName, title, leftLabel, rightLabel, cohortCompare, legendUpperRight)
{
  if (cohortCompare == F)     {     colorGrob1 <- "black"; colorGrob2 <- "black"     }
  if (cohortCompare == T)   { colorGrob1 <- "#ffa527"; colorGrob2 <- "#716be4"}
  annotationInfo <- paste0("NES: ", round(Pre$NES[grep(pathwayName,Pre$NAME)],2), "\n", "FDR: ", formatC(Pre$`FDR q-val`[grep(pathwayName,Pre$NAME)], format = "e", digits = 1))
  annotationInfo2 <- paste0("NES: ", round(Post$NES[grep(pathwayName,Post$NAME)],2), "\n", "FDR: ", formatC(Post$`FDR q-val`[grep(pathwayName,Post$NAME)], format = "e", digits = 1))
  if (legendUpperRight == T)   {   grob1 = grobTree(textGrob(annotationInfo, x=0.05,  y=0.52, hjust=0, gp=gpar(col="grey53", fontsize=20)))
  grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.22, hjust=0, gp=gpar(col="black", fontsize=20))) }
  if (legendUpperRight == F)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.25, hjust=0, gp=gpar(col="black", fontsize=24)))     }
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.97,hjust=0, gp=gpar(col=colorGrob1, fontsize=20)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.7,y=0.97,hjust=0, gp=gpar(col=colorGrob2, fontsize=20))) 
  return(
    ggplot(data=prePathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES) ) + geom_line(color="grey53", size=1) + 
      geom_rug(sides="t", size=0.75, alpha=0.5, color="grey53" ) + theme_bw() +
      geom_line(data=postPathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES),color="black", size=1)+ 
      geom_rug(data=postPathway, sides="b", size=0.75, alpha=0.5,color="black") +
      ggtitle(title) + ylab("Enrichment score") + xlab("Rank in gene list") + 
      theme(axis.text = element_text(size=19,hjust = 0.8, color="black"), axis.title = element_text(size=24,hjust = 0.5, color="black"), plot.title = element_text(size=32,hjust = 0.5))+
      annotation_custom(grob1) + annotation_custom(grob2)+ geom_hline(yintercept = 0) + 
      annotation_custom(left_grob) + annotation_custom(right_grob) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
      geom_hline(yintercept = 0, size=0.5)
  )
}

MitoticS_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_MITOTIC_SPINDLE.tsv", sep="\t")
MitoticS_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_MITOTIC_SPINDLE.tsv", sep="\t")

c<- plotGSEAvaccine(GSEA_PreBoost, GSEA_PostBoost, MitoticS_pre,MitoticS_post, pathwayName = "Mitotic_spindle", title = "Mitotic Spindle", 
                    leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
c

G2M_pre <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq/HALLMARK_G2M_CHECKPOINT.tsv", sep="\t")
G2M_post <- read.csv("~/S3data/GSEA_analyses/GSEA_analyses/Experienced_UP_Naive_DOWN_CD4_PostBoost_PseudoBulk_DESeq/HALLMARK_G2M_CHECKPOINT.tsv", sep="\t")

d<- plotGSEAvaccine(GSEA_PreBoost, GSEA_PostBoost, G2M_pre,G2M_post, pathwayName = "G2m_checkpoint", title = "G2M checkpoint", 
                    leftLabel = "", rightLabel = "", cohortCompare =F,legendUpperRight = T)
d

ggsave(a, filename = "./Images/Images/GSEA_IFNa_PostBooster.pdf", width = 7.25, height = 4.5)
ggsave(b, filename = "./Images/Images/GSEA_IFNg_PostBooster.pdf", width = 7.25, height = 4.5)
ggsave(c, filename = "./Images/Images/GSEA_MitoticS_PostBooster.pdf", width = 7.25, height = 4.5)
ggsave(d, filename = "./Images/Images/GSEA_G2M_PostBooster.pdf", width = 7.25, height = 4.5)

#--------pseudo bulk of AIM reactive CD4 subset on just pre booster timepoint--------------------
DefaultAssay(AIM.reactive.CD4) <- "RNA"
AIM.reactive.CD4.preBoost <- subset(AIM.reactive.CD4, timepoint == "7moPostVax" )
DefaultAssay(AIM.reactive.CD4.preBoost) <- "RNA"

#AIM.reactive.CD4.pre_post_boost <- subset(AIM.reactive.CD4, timepoint !=  "PostBreakThru" )
#View(AIM.reactive.CD4.pre_post_boost@meta.data)
AIM.reactive.CD4.sce <- AIM.reactive.CD4.preBoost %>% as.SingleCellExperiment() #using RNA counts !!!
View(AIM.reactive.CD4.sce)
View(AIM.reactive.CD4)
## Determine the number of cells per sample
table(AIM.reactive.CD4.sce$subject_timepoint)
groups <- colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")]
AIM.reactive.CD4.sce
AIM.reactive.CD4.sce <- removeAltExps(AIM.reactive.CD4.sce) 
# Aggregate across cluster-sample groups
pseudo_bulk_CD4 <- scuttle::aggregateAcrossCells(AIM.reactive.CD4.sce, ids = colData(AIM.reactive.CD4.sce)[, c("subject_timepoint")])

class(pseudo_bulk_CD4)
dim(pseudo_bulk_CD4)

#create metaData file 
metaData <- colnames(pseudo_bulk_CD4) %>% as.data.frame()
colnames(metaData) <- "subject_timepoint"
metaData <- metaData %>%
  separate(subject_timepoint,sep= "_",remove=FALSE, into=c("subject","timepoint"))
metaData <- metaData %>%
  mutate(PriorCOVID = case_when(
    endsWith(subject, "Subject1")  ~ "no",endsWith(subject, "Subject2")  ~ "no", endsWith(subject, "Subject3")  ~ "no",
    endsWith(subject, "Subject4")  ~ "no", endsWith(subject, "Subject5")  ~ "no",endsWith(subject, "Subject6")  ~ "no",
    endsWith(subject, "Subject7")  ~ "no",endsWith(subject, "Subject8")  ~ "yes",endsWith(subject, "Subject9") ~ "yes",
    endsWith(subject, "Subject10")  ~ "yes",endsWith(subject, "Subject11")  ~ "yes",endsWith(subject, "Subject12") ~ "yes",
    endsWith(subject, "Subject13") ~ "yes", endsWith(subject, "Subject14") ~ "yes", endsWith(timepoint, "Thru")  ~ "BreakThru"))
metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject, "Subject2")  ~ "NOV",endsWith(subject, "Subject3")  ~ "JAN",
    endsWith(subject, "Subject5")  ~ "JAN", endsWith(subject, "Subject6")  ~ "JAN",
    endsWith(subject, "Subject7")  ~ "DEC",endsWith(subject, "Subject8")  ~ "MARCH",
    endsWith(subject, "Subject9") ~ "DEC", endsWith(subject, "Subject12") ~ "MARCH",
    endsWith(subject, "Subject13") ~ "JAN",endsWith(subject, "Subject14") ~ "MARCH",
    endsWith(subject, "Subject1")  ~ "SEPT",endsWith(subject, "Subject4")  ~ "SEPT", endsWith(subject, "Subject5")  ~ "SEPT",
    endsWith(subject, "Subject10")  ~ "SEPT",endsWith(subject, "Subject11")  ~ "SEPT"))
metaData <- metaData %>%
  unite('subgroup', timepoint, PriorCOVID, remove=FALSE)
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
                              design = ~ subgroup) 
#subgroup "timepoint_PriorCOVID"

#Variance stabilizing transformation 
vsd <- vst(dds, blind=FALSE)
vsd$subgroup <- factor(vsd$subgroup, levels = c("7moPostVax_no","PostBreakThru_no","7moPostVax_yes"))
z <- plotPCA(vsd, intgroup = "runDate")
z + geom_point( size = 6)+ theme_classic()  

z <- plotPCA(vsd, intgroup = "subgroup")
z + geom_point( size = 6)+ theme_classic() + scale_color_manual(values = c("#ffc471", "#b5b2f1" )) 

#Batch variation removed using removeBatchEffect -- 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~subgroup, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design =mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "subgroup")

z <- plotPCA(vsd, intgroup = "runDate")
z + geom_point( size = 6)+ theme_classic()  

z <- plotPCA(vsd, intgroup = "subgroup")
z + geom_point( size = 2)+ theme_classic() + scale_color_manual(values = c("#ffc471","#b5b2f1" )) 


vst_df <- assay(vsd) %>%
  as.matrix()# Make into a data frame
#View(vst_df)
hist(vst_df)

View(YesvsNo_PREBoost)
#heatmap showing loadings 
#View(pc_loadings)
head(YesvsNo_PREBoost)
YesvsNo_PREBoost<- rownames_to_column(YesvsNo_PREBoost, var = "genes")
top.yes <-  YesvsNo_PREBoost %>% arrange((stat)) %>% head(15)
top.no <-  YesvsNo_PREBoost %>% arrange(-(stat)) %>% head(15)
top <- c(top.yes$genes,top.no$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() #22,133 genes of 22 samples
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "PriorCOVID")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  PriorCOVID = c(yes = "#b5b2f1", no = "#ffc471" ))
z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = T,
                        scale = 'row', cutree_cols = 1,treeheight_row = FALSE,treeheight_col = TRUE,
                        border_color = NA,
                        annotation_col = metaData, # Add metadata labels!
                        show_colnames = T, # Don't show sample labels
                        fontsize = 10, # Shrink the pathway labels
                        annotation_colors	= cols, 
                        color = colorRampPalette(c("blue", "white", "firebrick3"))(50))

z
ggsave(z, filename = "./Images/Images/Heatmap_PriorCOVID_rankedDEseq_PreBoost.pdf", width = 4, height = 5)
ggsave(z, filename = "./Images/Images//Heatmap_PriorCOVID_rankedDEseq_PreBoost.jpeg", width = 4, height = 5)

pseudo.bulk <- assay(vsd) %>% as.data.frame() #
head(pseudo.bulk)
pseudo.bulk <- pseudo.bulk %>% as.matrix() %>% t()

pseudo_bulk_pca <- prcomp(pseudo.bulk)
pc_eigenvalues <- pseudo_bulk_pca$sdev^2

pc_scores <- pseudo_bulk_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores
pc_scores <- merge(pc_scores, metaData, by = 1)
pc_scores$PriorCOVID <- factor(pc_scores$PriorCOVID, levels = c("no","yes"))
pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, color = PriorCOVID)) +
  geom_point(size = 4)+ theme_classic() + scale_color_manual(values = c("#ffc471","#b5b2f1"))

pc_loadings <- pseudo_bulk_pca$rotation
pc_loadings <- pc_loadings %>% 
  as.data.frame()
pc_loadings$genes <- row.names(pc_loadings)
# print the result
head(pc_loadings)
#View(pc_loadings)
row.names(pc_loadings)
top_PC1 <-  pc_loadings %>% arrange(-abs(PC1)) %>% head(20)
top_PC2 <-  pc_loadings %>% arrange(-abs(PC2)) %>% head(20) 

top <- c(top_PC1$genes, top_PC2$genes)
pseudo.bulk <- assay(vsd) %>% as.data.frame() 
head(pseudo.bulk)
pseudo.bulk.top <- pseudo.bulk %>% 
  filter(rownames(pseudo.bulk) %in% top) %>% 
  as.matrix

metaData <- metaData[,c("subject_timepoint", "PriorCOVID")]
metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("subject_timepoint")
cols <-  list(
  PriorCOVID = c(yes = "#b5b2f1", no = "#ffc471" ))

z <- pheatmap::pheatmap(pseudo.bulk.top, cluster_rows=T,cluster_cols = TRUE,
                   scale = 'row',cutree_cols = 1,treeheight_row = FALSE,treeheight_col = FALSE,
                   border_color = NA,
                   annotation_col = metaData, # Add metadata labels!
                   show_colnames = FALSE, # Don't show sample labels
                   fontsize_row = 13, # Shrink the pathway labels
                   annotation_colors	= cols, fontsize = 13,
                   color = colorRampPalette(c("blue", "white", "firebrick3"))(50))  
z 

ggbiplot(pseudo_bulk_pca, ellipse = TRUE, groups = metaData$PriorCOVID,                   
         var.axes = FALSE )  +  # Remove variable vectors (TRUE)
        theme_classic() + scale_colour_manual(name = "Prior COVID?",values= c( "#ffc471", "#b5b2f1"))+
  labs(title="pseudobulk of AIM reactive CD4+ T cells (pre boost)")

#DE Analysis
DESdata_mixed <- DESeq(dds, normalized=FALSE, parallel=TRUE)

data <- as.data.frame(results(DESdata_mixed, contrast=c("PriorCOVID","yes","no"), parallel = T))
View(data)
data$diffexpressed <- "NO"
data$diffexpressed[data$log2FoldChange > 0.2 & 
                     data$pvalue <= 0.05] <- "UP"
data$diffexpressed[data$log2FoldChange < -0.2 & 
                     data$pvalue <= 0.05] <- "DOWN"
UP <- subset(data, subset = diffexpressed == "UP" )
DOWN <- subset(data, subset = diffexpressed == "DOWN" )
top <- rbind(UP,DOWN)%>% row.names()
top
plot2<- ggplot(data=data, aes(x=log2FoldChange, y=-log10(pvalue), col =diffexpressed))+
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "lightgrey")+
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "longdash", color="lightgrey")+
  geom_point()+ theme_classic()+ 
  scale_color_manual(values = c("#ffc471","lightgrey","#b5b2f1"))
z <-LabelPoints(plot = plot2, points = top, repel = TRUE, xnudge = 0, ynudge = 0, size = 4)+labs(title="AIM Reactive CD4+ T cells (pseudobulk, pre boost)")+
  NoLegend()+labs(y = "-log10 p value", x = "Naive -------------- log2 (Fold Change) -------------- Experienced")+
  theme(panel.background=element_rect(fill="white"), 
        panel.border=element_rect(fill = NA, colour = "black", size = 1))
z
write.csv(data,"./Gene_sets/Experienced_UP_Naive_DOWN_CD4_PreBoost_PseudoBulk_DESeq.csv", row.names = TRUE)

#-------------------------------MAPPING TO ADAPTIVE DATABASE-------------------------------
#split CTaa into CT_TRA and CT_TRB
View(AIM.integrated.T.clean.CD4@meta.data)
table(AIM.integrated.T.clean.CD4$integrated_cluster.IDs)
AIM.integrated.T.clean.CD4@meta.data <- separate(AIM.integrated.T.clean.CD4@meta.data, 
                                       CTaa, sep= "_",remove=FALSE, into=c("CT_TRA","CT_TRB"))
AIM.integrated.T.clean.CD4@meta.data$CT_TRB[AIM.integrated.T.clean.CD4@meta.data$CT_TRB == "NA"] <- NA

AIM.integrated.T.clean.CD4@meta.data$TCR_B <- !is.na(AIM.integrated.T.clean.CD4@meta.data$CT_TRB)
View(AIM.integrated.T.clean.CD4@meta.data)

DimPlot(AIM.integrated.T.clean.CD4)
AIM.integrated.T.clean.CD4.TCRB <- subset(AIM.integrated.T.clean.CD4, TCR_B == "TRUE")

z<- DimPlot(AIM.integrated.T.clean.CD4, group.by = "cloneType", raster = F, order = c("Hyperexpanded (0.01 < X <= 1)",
                                                                                      "Large (0.005 < X <= 0.01)",
                                                                                      "Medium (0.001 < X <= 0.005)",
                                                                                      "Small (5e-04 < X <= 0.001)",
                                                                                      "Rare (0 < X <= 5e-04)")) +
  scale_color_manual(values = colorblind_vector(5), na.value=NA) + 
  theme(plot.title = element_blank())
z
ggsave(z, filename = "./Images/Images/CD4_UMAP_CloneType.jpeg", width = 8, height = 4)
ggsave(z, filename = "./Images/Images/CD4_UMAP_CloneType.pdf", width = 8, height = 4)

#Class II 
peptide.detail.cii <- read_csv("~/S3data/Adaptive_peptide_detail/peptide-detail-cii.csv") #6809 rows 
peptide.detail.cii$Class <- "2"
peptide.detail.cii <- peptide.detail.cii %>% separate(`TCR BioIdentity`, c("cdr3", "v_gene", "j_gene"), "\\+")
peptide.detail.cii.surface <- peptide.detail.cii[grep("surface glycoprotein", peptide.detail.cii$`ORF Coverage`), ]
length(unique(peptide.detail.cii.surface$cdr3))
View(peptide.detail.cii.surface)
AIM.integrated.T.clean.CD4.TCRB@meta.data <-  AIM.integrated.T.clean.CD4.TCRB@meta.data %>% 
  add_column(MapsToAdaptive = 
               AIM.integrated.T.clean.CD4.TCRB$CT_TRB %in% peptide.detail.cii.surface$cdr3)

z<- DimPlot(AIM.integrated.T.clean.CD4.TCRB, group.by = "MapsToAdaptive", pt.size = 1, order = c("TRUE","FALSE"), 
            split.by = "integrated_cluster.IDs",
        cols = c("lightgray","#F68282")) +labs(title = "CD4+ T cdr3B sequences that map to adaptive")
z
ggplot(AIM.integrated.T.clean.CD4.TCRB@meta.data, aes(fill=factor(MapsToAdaptive), x=integrated_cluster.IDs,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") + 
  scale_fill_manual(values = c("lightgray","#F68282"))

AIM.integrated.T.clean.CD4.TCRB$integrated_cluster.IDs <- factor(AIM.integrated.T.clean.CD4.TCRB$integrated_cluster.IDs,
                                                levels = c("AIM Reactive CD4+ T","CD4+ T","Naive CD4+ T","Treg"))
data <- table(AIM.integrated.T.clean.CD4.TCRB$integrated_cluster.IDs, AIM.integrated.T.clean.CD4.TCRB$subject, 
              AIM.integrated.T.clean.CD4.TCRB$MapsToAdaptive) %>% as.data.frame()
colnames(data) <- c("cluster","subject","overlaps","count")
false <- subset(data, overlaps == "FALSE")
true <- subset(data, overlaps == "TRUE")
data <- merge(false,true, by = c("cluster","subject"))
data$percentOverlap <- data$count.y/sum(data$count.x,data$count.y)*100
data$uniqueAdaptiveMINUSoverlap <- 2783-data$count.y
data <- data %>% unite("subject_cluster",subject,cluster,remove=F)
data.list <- setNames(split(data, seq(nrow(data))), data$subject_cluster)
View(data.list)
for (i in seq_along(data.list)) {
  data.list[[i]] <- data.frame(repeating=c(data.list[[i]]$count.y,(data.list[[i]]$uniqueAdaptiveMINUSoverlap)),
                               notRepeat=c(data.list[[i]]$count.x,10000000))
}
for (i in seq_along(data.list)) {
  data.list[[i]]  <- -log10(fisher.test(data.list[[i]])$p.value)
}
data.list.df <- as.data.frame(do.call(rbind, data.list))
colnames(data.list.df) <- "logpvalueFisher"
data.list.df <- rownames_to_column(data.list.df, "subject_cluster")

data.2 <- merge(data,data.list.df, by = "subject_cluster")
View(data.2)
data.2$cluster <- factor(data.2$cluster, levels = c("AIM Reactive CD4+ T","CD4+ T","Naive CD4+ T","Treg"))
levels(data.2$cluster)
z<-ggplot(data.2,aes(x=cluster,y=percentOverlap, group = subject)) +geom_jitter(aes(size=logpvalueFisher),width = 0.25, height = 0)+ 
   scale_size(range = c(0.5, 5)) + 
  theme(legend.key = element_rect(fill=NA)) + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y=element_line(colour="gray",linetype="dashed"),
        axis.text=element_text(colour="black",size=15))  + ylim(-0.001,0.025)+labs(title = "CD4 overlap with Adaptive Class ii cdr3") +
  stat_summary(aes(group = cluster), fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5, size = 0.25)
z
ggsave(z, filename="./Images/Images/clonalOverlap_PerCluster_CD4_classii_Adaptive.conservative.per.subject.pdf", width = 5, height = 5)

data <- table(AIM.integrated.T.clean.CD4.TCRB$integrated_cluster.IDs, 
              AIM.integrated.T.clean.CD4.TCRB$MapsToAdaptive) %>% as.data.frame()
colnames(data) <- c("cluster","overlaps","count")
false <- subset(data, overlaps == "FALSE")
true <- subset(data, overlaps == "TRUE")
data <- merge(false,true, by = c("cluster"))
data$percentOverlap <- data$count.y/sum(data$count.x,data$count.y)*100
data$uniqueAdaptiveMINUSoverlap <- 2783-data$count.y
data.list <- setNames(split(data, seq(nrow(data))), data$cluster)
View(data.list)
for (i in seq_along(data.list)) {
  data.list[[i]] <- data.frame(repeating=c(data.list[[i]]$count.y,(data.list[[i]]$uniqueAdaptiveMINUSoverlap)),
                               notRepeat=c(data.list[[i]]$count.x,10000000))
}
for (i in seq_along(data.list)) {
  data.list[[i]]  <- -log10(fisher.test(data.list[[i]])$p.value)
}
data.list.df <- as.data.frame(do.call(rbind, data.list))
colnames(data.list.df) <- "logpvalueFisher"
data.list.df <- rownames_to_column(data.list.df, "cluster")

data.2 <- merge(data,data.list.df, by = "cluster")
View(data.2)
data.2$cluster <- factor(data.2$cluster, levels = c("AIM Reactive CD4+ T","CD4+ T","Naive CD4+ T","Treg"))
levels(data.2$cluster)
z<-ggplot(data.2,aes(x=cluster,y=percentOverlap, group = cluster)) +geom_jitter(aes(size=logpvalueFisher),width = 0, height = 0)+ 
  scale_size(range = c(1, 10)) + 
  theme(legend.key = element_rect(fill=NA)) + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y=element_line(colour="gray",linetype="dashed"),
        axis.text=element_text(colour="black",size=12))  + ylim(0,0.15)+labs(title = "CD4 overlap with Adaptive Class ii cdr3")
z
ggsave(z, filename="./Images/Images/clonalOverlap_PerCluster_CD4_classii_Adaptive.conservative.pdf", width = 5, height = 5)

#Just AIM reactive clusters mapping to Adaptive and graphed per subject
#AIM reactive CD4s
AIM.CD4 <- subset(AIM.integrated.T.clean.CD4.TCRB, integrated_cluster.IDs == "AIM Reactive CD4+ T")
table(AIM.CD4$MapsToAdaptive) 
ggplot(AIM.CD4@meta.data, aes(fill=factor(MapsToAdaptive), x=subject,)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "freq", x = "Cluster") + 
  scale_fill_manual(values = c("lightgray","#F68282"))
length(unique(peptide.detail.cii.surface$cdr3))
length(unique(AIM.CD4$CT_TRB))

data <- table(AIM.CD4$subject, 
              AIM.CD4$MapsToAdaptive) %>% as.data.frame()
colnames(data) <- c("subject","overlaps","count")
false <- subset(data, overlaps == "FALSE")
true <- subset(data, overlaps == "TRUE")
data <- merge(false,true, by = c("subject"))
data$percentOverlap <- data$count.y/sum(data$count.x,data$count.y)*100
data$uniqueAdaptiveMINUSoverlap <- 2783-data$count.y
data.list <- setNames(split(data, seq(nrow(data))), data$subject)

for (i in seq_along(data.list)) {
  data.list[[i]] <- data.frame(repeating=c(data.list[[i]]$count.y,(data.list[[i]]$uniqueAdaptiveMINUSoverlap)),
                               notRepeat=c(data.list[[i]]$count.x,10000000))
}
for (i in seq_along(data.list)) {
  data.list[[i]]  <- -log10(fisher.test(data.list[[i]])$p.value)
}
data.list.df <- as.data.frame(do.call(rbind, data.list))
colnames(data.list.df) <- "logpvalueFisher"
data.list.df <- rownames_to_column(data.list.df, "subject")

data.2 <- merge(data,data.list.df, by = "subject")
View(data.2)
data.2$PriorCOVID <- "yes"
data.2[1:7,"PriorCOVID"] <- "no"
data.2$subject <- factor(data.2$subject , levels = data.2$subject [order(data.2$percentOverlap, decreasing=TRUE)])

z<-ggplot(data.2,aes(x=subject,y=percentOverlap, group = subject, color = PriorCOVID)) +
  geom_jitter(aes(size=logpvalueFisher),width = 0, height = 0)+ scale_color_manual(values = c(vax, infect))+
  scale_size(range = c(1, 7)) + 
  theme(legend.key = element_rect(fill=NA), legend.text = element_text(size = 15)) + theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y=element_line(colour="gray",linetype="dashed"),
        axis.text=element_text(colour="black",size=10))  + ylim(0,0.2)+labs(title = "CD4 overlap with Adaptive Class ii cdr3")
z
ggsave(z, filename = "./Images/Images/clonalOverlap_AIM_CD4_classii_Adaptive.conservative.subject.pdf", width =5 , height = 3)

#looking at antigens recognized
z <- merge(AIM.CD4@meta.data, peptide.detail.cii.surface, by.x = "CT_TRB", by.y = "cdr3")
S_protein_positions <- read_csv("Adaptive_peptide_detail/S_protein_positions.csv")
x <- merge(z, S_protein_positions, by.x = "Start Index in Genome", by.y = "...1")
x <- merge(x, S_protein_positions, by.x = "End Index in Genome", by.y = "...1")
x <- x %>% unite("S_Protein_Postion", `...2.x`, `...2.y`, sep = "-")
View(x)
overlap <- x[,c("Start Index in Genome", "End Index in Genome", "CT_TRB", "S_Protein_Postion")]
ggplot(overlap, aes(fill=factor(S_Protein_Postion), x="",)) + 
  geom_bar(position="fill")+theme_classic()+labs(y = "freq") 

freq <- prop.table(table(overlap$S_Protein_Postion)) %>% as.data.frame()
freq <- freq %>% arrange(-Freq)
freq$Freq <- signif(freq$Freq, 3)
library(ggrepel)
freq.sub <- subset(freq, Freq > 0.05)
freq2 <- freq.sub %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

coul <- brewer.pal(9, name = "Reds") 
coul <- colorRampPalette(coul)(13)

z <- ggplot(freq, aes(x="", y=Freq, fill=reorder(Var1, Freq))) +
  geom_bar(stat="identity",color="white") +
  coord_polar("y") + coord_polar(theta = "y") +
  geom_label_repel(data = freq2,
                   aes(y = pos, label = paste0(Freq*100, "%")),
                   size = 3, nudge_x = 0.75, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Spike peptide position")) +
  theme_void()+ scale_fill_manual(values = rev(coul))
z
ggsave(z, filename = "./Images/Images//Adaptive_overlap_peptide_positions.pdf", width =5 , height = 4)


#---------------------Highlight clones------------------
Idents(AIM.integrated.T.clean.CD4) <- "integrated_cluster.IDs"
combined.AIM.CD4 <- expression2List(AIM.integrated.T.clean.CD4, split.by = "integrated_cluster.IDs")
x<- compareClonotypes(combined.AIM.CD4,number = 21, samples = c("AIM Reactive CD4+ T"), 
                      cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
View(x) #look at top 20 clones 
top <- x %>% head(20)
top_20_clones <- top$Clonotypes 
top_20_clones <- as.character(top_20_clones)
top_20_clones

seurat <- highlightClonotypes(AIM.integrated.T.clean.CD4, 
                              cloneCall= "aa",
                              sequence = top_20_clones)
                                
seurat@meta.data$highlight <- as.factor(ifelse(is.na(seurat@meta.data$highlight), "NA", seurat@meta.data$highlight))
View(seurat@meta.data)

colorblind_vector <- colorRampPalette((c("lightgrey","#0D0887FF", "#0D0887FF", "#7301A8FF", 
                                                    "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                                    "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
                                                    

z <- DimPlot(seurat, group.by = "highlight", raster = FALSE,pt.size = 0.5, order = c(top_20_clones, "NA")) + 
  theme(plot.title = element_blank())+  scale_color_manual(values = colorblind_vector(21),na.value = "lightgrey") 
z


ggsave(z,filename = "~/S3data/Images/Images/Top10Clones_AIM_Reactive_cluster.jpeg",width =13,height = 4)
ggsave(z,filename = "~/S3data/Images/Images/Top10Clones_AIM_Reactive_cluster.pdf",width = 13,height = 4)

#-----------GSVA infection and vaccine imprints------------
AIM.reactive.CD4@meta.data <- AIM.reactive.CD4@meta.data %>% unite('COVID_timepoint',PriorCOVID, timepoint,  remove=FALSE)
DefaultAssay(AIM.reactive.CD4) <- "RNA"
#AIM.reactive.CD4.preBoost <- NormalizeData(AIM.reactive.CD4.preBoost)
ActivatedT.mxt.RNA <- as.data.frame(AIM.reactive.CD4@assays$RNA@scale.data) #36601 genes
#View(ActivatedT.mxt)
ActivatedT.mxt.RNA <- as.matrix(ActivatedT.mxt.RNA)
#View(ActivatedT.mxt.RNA)

imprints <- list(rownames(Infection_imprint),rownames(Vaccine_imprint))
names(imprints) = c("Infection_imprint","Vaccine_imprint")
View(imprints)
gsva.es.all <- gsva(ActivatedT.mxt.RNA, imprints, verbose=TRUE, method="gsva")

#saveRDS(gsva.es.all, file = "~/S3data/Saved_Objects/gsva_infection_vax_imprints.RDS") 
gsva.es.all <- readRDS("~/S3data/Saved_Objects/gsva_infection_vax_imprints.RDS")

dim(gsva.es.all)

metaData <- colnames(ActivatedT.mxt.RNA)
data_to_add <- rownames_to_column(AIM.reactive.CD4@meta.data ,"barcodes")
data_to_add <- data_to_add[,c("barcodes","subject","timepoint","PriorCOVID", "COVID_timepoint")]
colnames(data_to_add) <- c("barcode","subject","timepoint","PriorCOVID", "COVID_timepoint") 

metaData <- merge(metaData, data_to_add, by.x = 1, by.y = 1)
## Determine the number of cells per sample

colnames(metaData) <- c("barcode", "subject","timepoint","PriorCOVID", "COVID_timepoint") 
metaData <- metaData %>% arrange(subject)

gsva.es.all <- gsva.es.all[ , metaData$barcode]


metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("barcode")
cols <-  list(
  PriorCOVID = c(yes = "#b5b2f1", no = "#ffc471" ))

pathway_heatmap <- pheatmap::pheatmap(gsva.es.all,annotation_col = metaData, # Add metadata labels!
                                      #cutree_cols = 2,
                                      scale = "row",
                                      cluster_cols = T,
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6, # Shrink the pathway labels
                                      #annotation_colors	= cols, 
                                      color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

gsva.results <- t(gsva.es.all)
data <- merge(metaData, gsva.results, by.x = 0, by.y = 0) %>% as.data.frame()

ggplot(data, aes(x=subject,y=Infection_imprint, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Infection imprint") + facet_grid(rows="timepoint")

ggplot(data, aes(x=subject,y=Vaccine_imprint, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Vaccine imprint") + facet_grid(rows="timepoint")

data$timepoint <- factor(data$timepoint, levels = c("7moPostVax","1moPostBoost","PostBreakThru"))
data$COVID_timepoint <- factor(data$COVID_timepoint, levels = c("no_7moPostVax","no_1moPostBoost","no_PostBreakThru","yes_7moPostVax","yes_1moPostBoost"))
z <- ggplot(data, aes(x=Infection_imprint,y=Vaccine_imprint, color = COVID_timepoint)) +
  theme_classic() + geom_point(size=0.05) + geom_density2d(binwidth = 0.25,color="grey", size = 0.25) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(pre_vax,post_vax,breakthru,pre_infect,post_infect)) + labs(x =  "Infection signature enrichment score",
                                                                                           y="Vaccine signature enrichment score") +  
  facet_grid(PriorCOVID ~ timepoint) + NoLegend() + theme(line=element_line(colour = "black", size = 0.25)) +
 scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots.pdf",width = 6,height = 4)
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots.jpeg",width = 6,height = 4)

View(data)

z <- ggplot(data, aes(x=Infection_imprint,y=Vaccine_imprint, color = COVID_timepoint)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,  size = 0.15)+ geom_point(size=0.05) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(pre_vax,post_vax,breakthru,pre_infect,post_infect)) + labs(x =  "Infection signature enrichment score",
                                                                                           y="Vaccine signature enrichment score") +  
  facet_grid(subject ~ timepoint) + NoLegend()+ theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_individual.pdf",width = 8,height = 18)
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_individual.jpeg",width = 8,height = 18)

table(AIM.reactive.CD4$subject, AIM.reactive.CD4$timepoint)

View(test)
data <- data %>%
  mutate(Signature = case_when(
    Vaccine_imprint > 0 & Infection_imprint < 0   ~ "Vaccine_Signature",
    Vaccine_imprint < 0 & Infection_imprint > 0   ~ "Infection_Signature", 
    Vaccine_imprint > 0 & Infection_imprint > 0 ~ "Double Positive",
    Vaccine_imprint < 0 & Infection_imprint < 0 ~ "Double Negative"))
vaccine_pos <- subset(data, Vaccine_imprint > 0)

infection_pos <- subset(data, Infection_imprint > 0)
data$Signature <- factor(data$Signature, levels = c("Double Positive","Vaccine_Signature","Infection_Signature","Double Negative"))
z <- ggplot(data, aes(fill=factor(Signature), x=timepoint)) +
  theme_classic() +  geom_bar(position="fill") +
  scale_fill_manual(values = c( "lightgrey",vax,infect,"black")) + labs(x =  "",
                                                                        y="Count") +
  facet_grid(~PriorCOVID)
z 

#frequencies in each qunadrants
Signature_groups <- data %>% 
  group_by(subject,timepoint,Signature) %>%
  dplyr::summarize(count= n()) %>% dplyr::ungroup() %>% as.data.frame()

Total <- data %>% 
  group_by(subject,timepoint) %>%
  dplyr::summarize(total= n()) %>% dplyr::ungroup() %>% as.data.frame()

df <- merge(Signature_groups,Total, by = 1:2)
df$freq <- df$count/df$total
df.wide <- df[,c("subject","timepoint","Signature","freq")]
df.wide <- spread(df.wide, Signature, freq)
View(df.wide)
#write.csv(df,"./Export_csv/ContourPlotsFrequencies.csv", row.names = TRUE)
#write.csv(df.wide,"./Export_csv/ContourPlotsFrequencies.wide.csv", row.names = TRUE)

#sum of vaccine and infection signatures 
Infection_score <- data %>% 
  group_by(subject,timepoint, Signature) %>%
  dplyr::summarize(sum(Infection_imprint)) %>% dplyr::ungroup() %>% as.data.frame()

Vaccine_score <- data %>% 
  group_by(subject,timepoint,Signature) %>%
  dplyr::summarize(sum(Vaccine_imprint)) %>% dplyr::ungroup() %>% as.data.frame()

df <- merge(Infection_score,Vaccine_score, by = 1:3)
View(df)
df$`sum(Infection_imprint)`

ggplot(df, aes(x=subject,y=`sum(Vaccine_imprint)`)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Infection imprint") + facet_grid(rows="Signature")


gsva.results.1 <- gsva.results %>% as.data.frame() %>% rownames_to_column( "barcodes")
AIM.reactive.CD4@meta.data  <- rownames_to_column(AIM.reactive.CD4@meta.data, "barcodes")
AIM.reactive.CD4@meta.data <- merge(AIM.reactive.CD4@meta.data, gsva.results.1, by = "barcodes")
AIM.reactive.CD4@meta.data <- column_to_rownames(AIM.reactive.CD4@meta.data, "barcodes")
View(AIM.reactive.CD4@meta.data)

AIM.reactive.CD4@meta.data$timepoint <- factor(AIM.reactive.CD4@meta.data$timepoint , levels = c("7moPostVax" , "1moPostBoost", "PostBreakThru"))
z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=Infection_imprint,y=Vaccine_imprint, color = RNA_snn_res.0.4)) +
  theme_classic()  + geom_density2d(binwidth = 0.25,color="lightgrey", size = 0.25) +geom_point(size=0.25)+
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_brewer(palette = "Paired") + labs(x =  "Infection signature enrichment score",
                                                                y="Vaccine signature enrichment score") +  
  facet_grid(RNA_snn_res.0.4 ~ COVID_timepoint) + theme(line=element_line(colour = "black", size = 0.25))

z
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType.pdf",width = 6,height = 4)
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType.jpeg",width = 6,height = 4)


combined.AIM.CD4 <- expression2List(AIM.reactive.CD4, split.by = "subject_timepoint")
x<- compareClonotypes(combined.AIM.CD4,number = 5,  
                      cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
View(x) #look at top 20 clones 
top <- x #%>% head(20)
top_20_clones <- top$Clonotypes 
top_20_clones <- as.character(top_20_clones)
top_20_clones

seurat <- highlightClonotypes(AIM.reactive.CD4, 
                              cloneCall= "aa",
                              sequence = c("CAVRYSLNSGYSTLTF_CASSELLAGGPNEQYF"))
#View(seurat@meta.data)
z <- ggplot(seurat@meta.data, aes(x=Infection_imprint,y=Vaccine_imprint, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="grey", size = 0.25)+ geom_point(size=3) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = rev(colorblind_vector(55)), na.value = NA) + labs(x =  "Infection signature enrichment score",
                                                                y="Vaccine signature enrichment score") +  
  facet_grid(PriorCOVID ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25)) + NoLegend()

z

AIM.reactive.CD4 <- highlightClonotypes(AIM.reactive.CD4, 
                              cloneCall= "aa", 
                              sequence = c(#"CAADFSGGGADGLTF_CASRTGSTNEKLFF", 
                                           #"CAMRAYNTGNQFYF_CASTSQGASDTEAFF",
                                          #"CAMRPLNTGNQFYF_CASSQESAGGIDEQFF",
                                         #  "CAVYTSGTYKYIF_CSVDGQGNTGELFF", 
                                          "CAVQDLLASGSRLTF_CASSSPSGELFF","CALRDGGSNYKLTF_CAWSASPDTQYF", #Subject1
                                          "CALLANQAGTALIF_CASSESNRGAYGYTF","CAVQALLNAGNMLTF_CASGARTGEQFF", #Subject1
                                          "CAVVPLGNTGKLIF_CASSQDLTGANVLTF", #Subject4
                                          "CAGAEETSGSRLTF_CASSLSGGSNSPLHF", #Subject5
                                          "CAASEAAGNKLTF_CASSQVGGMGAKNIQYF", #Subject6
                                          #"CAMREVNTGTASKLTF_CASIPPGSAQGAFQPQHF",
                                          # "CAGRRGKLIF_CASSYEGPFGEQFF" ,
                                          "CAMRLGGSYIPTF_CASSPPRLQDNEQFF" , # Subject14
                                           "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF", # Subject14
                                          "CATDWNNDMRF_CASSFTGTGNEKLFF",   #Subject8 
                                          "CALNTGGFKTIF_CASSEAGGAGRVNEQFF"  #Subject11
                                                                                                                                                                           ))
AIM.reactive.CD4.vax <- subset(AIM.reactive.CD4, PriorCOVID == "no")
AIM.reactive.CD4.vax <- highlightClonotypes(AIM.reactive.CD4.vax, 
                                        cloneCall= "aa", 
                                        sequence = c("CALRDGGSNYKLTF_CAWSASPDTQYF",
                                          "CAVQDLLASGSRLTF_CASSSPSGELFF","CALLANQAGTALIF_CASSESNRGAYGYTF",#"CAVQALLNAGNMLTF_CASGARTGEQFF","CVVHGNNFNKFYF_CSARAGEGNYGYTF",#Subject1
                                          "CAVVPLGNTGKLIF_CASSQDLTGANVLTF" #Subject4 
                                         # "CAASQGAQKLVF_CASSQARLAPQFF",      #Subject5
                                         # "CAASEAAGNKLTF_CASSQVGGMGAKNIQYF", "CAAFISGNTPLVF_CSTGGRGGETQYF" #Subject7
                                          ))   
vax_clones = c( "#f6aa1c", "#f25c54", "#bf3100","#550527") #, "orange", "yellow2", "red", "red3","brown" )

z <- ggplot(AIM.reactive.CD4.vax@meta.data, aes(x=Infection_imprint,y=Vaccine_imprint, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="lightgrey", size = 0.25)+ geom_point(size=0.5) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values=vax_clones , na.value=NA) + labs(x =  "Infection signature enrichment score",
                                                                y="Vaccine signature enrichment score") +  
  facet_grid( ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25))+
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType_vax.pdf",width = 8,height = 2)
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType_vax.jpeg",width = 8,height = 2)

pre <- prop.table(table(AIM.reactive.CD4.vax$highlight, AIM.reactive.CD4.vax$RNA_snn_res.0.4), margin = 2)%>% as.data.frame()
pre$Freq <- pre$Freq*100
colnames(pre) <- c("Clone","Cluster", "%")
z <- ggplot(pre, aes(x=`Clone`, y=`%`, fill = Cluster, group = Cluster)) +
  geom_col(width = 0.75, colour = "white") +
  #geom_area(aes(x = c("no" = 1.25, "yes" = 1.75)[`PriorCOVID`]), 
  # position = "fill", colour = "white", alpha = 0.4,
  #outline.type = "both") + 
  theme_classic()  + scale_fill_brewer(palette = "Paired") +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="Prior COVID?", title = "Pre booster")+ RotatedAxis()
z



subset <- subset(AIM.reactive.CD4.vax, highlight == "CALRDGGSNYKLTF_CAWSASPDTQYF" & PriorCOVID == "no" |
                   highlight == "CAVQDLLASGSRLTF_CASSSPSGELFF"& PriorCOVID == "no" |highlight == "CALLANQAGTALIF_CASSESNRGAYGYTF" & PriorCOVID == "no"|
                   highlight == "CAVVPLGNTGKLIF_CASSQDLTGANVLTF" & PriorCOVID == "no" ) 

infection_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
    dplyr::summarize(sum(Infection_imprint)) 

vaccine_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Vaccine_imprint)) 

ggplot(vaccine_score_vax,aes(x=timepoint,y=`sum(Vaccine_imprint)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of vaccine signature enrichment score")
ggplot(infection_score_vax,aes(x=timepoint,y=`sum(Infection_imprint)`, color = highlight, group = highlight)) + geom_point()+geom_line()+ theme_classic()+ 
  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of infection signature enrichment score")

z <- ggplot(subset@meta.data, aes(x=timepoint,y=Vaccine_imprint, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
    geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+ 
  labs(x =  "cluster", y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight")  + scale_y_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Vaccine_imprints_VlnPlot_CloneType_vax.pdf",width = 5.25,height = 6)

z<- ggplot(subset@meta.data, aes(x=timepoint,y=Infection_imprint, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+
  labs(x =  "cluster", y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Infection_imprints_VlnPlot_CloneType_vax.pdf",width = 5.25,height = 6)

export <- subset@meta.data[,c("barcode","timepoint","Infection_imprint","highlight")]
export.wide <- export %>% unite("Clone_timepoint", highlight, timepoint, remove = T)
export.wide <- spread(export.wide, Clone_timepoint, Infection_imprint)
View(export.wide)
write.csv(export,"~/S3data/Export_csv/Vaccine_Top_Clones_GSVA_Scores.csv", row.names = F)

#infection-primed 
AIM.reactive.CD4.infect<- subset(AIM.reactive.CD4, PriorCOVID == "yes")
AIM.reactive.CD4.infect <- highlightClonotypes(AIM.reactive.CD4.infect, 
                                            cloneCall= "aa", 
                                            sequence = c(#"CATDWNNDMRF_CASSFTGTGNEKLFF", #Subject8
                                                         "CAMRLGGSYIPTF_CASSPPRLQDNEQFF",#Subject14, Subject12
                                                         "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF", #"Subject14
                                                         "CALNTGGFKTIF_CASSEAGGAGRVNEQFF", 
                                                         "CAADFSGGGADGLTF_CASRTGSTNEKLFF" #Subject11
                                                         #"CAVRYSLNSGYSTLTF_CASSELLAGGPNEQYF" #Subject14, Subject8  
                                                         ))
infect_clones = c( "#E53F71","#8338ec","#E090DF", "#541675")

z <- ggplot(AIM.reactive.CD4.infect@meta.data, aes(x=Infection_imprint,y=Vaccine_imprint, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="lightgrey", size = 0.25)+ geom_point(size=0.5) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = infect_clones, na.value=NA) + labs(x =  "Infection signature enrichment score",
                                                           y="Vaccine signature enrichment score") +  
  facet_grid( ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType_infect.pdf",width = 6.5,height = 2)
ggsave(z,filename = "~/S3data/Images/Images/Infection_Vaccine_imprints_ContourPlots_CloneType_infect.jpeg",width = 6.5,height = 2)

subset <- subset(AIM.reactive.CD4.infect, #highlight ==   "CATDWNNDMRF_CASSFTGTGNEKLFF" & PriorCOVID == "yes"|   #Subject8
                 highlight ==    "CAMRLGGSYIPTF_CASSPPRLQDNEQFF" & PriorCOVID == "yes"|  # Subject14
                 highlight ==   "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF" & PriorCOVID == "yes" | # Subject14
                 highlight ==    "CALNTGGFKTIF_CASSEAGGAGRVNEQFF" & PriorCOVID == "yes"|
                   highlight ==    "CAADFSGGGADGLTF_CASRTGSTNEKLFF" & PriorCOVID == "yes"
                  # highlight ==    "CAVRYSLNSGYSTLTF_CASSELLAGGPNEQYF" & PriorCOVID == "yes" 
                   )

infection_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Infection_imprint)) 

vaccine_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Vaccine_imprint)) 

ggplot(vaccine_score_vax,aes(x=timepoint,y=`sum(Vaccine_imprint)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values = infect_clones, na.value=NA) +labs(x="time point", y="Sum of Vaccine Imprint")
ggplot(infection_score_vax,aes(x=timepoint,y=`sum(Infection_imprint)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values = infect_clones, na.value=NA) +labs(x="time point", y="Sum of Infection Imprint")

z<-ggplot(subset@meta.data, aes(x=timepoint,y=Vaccine_imprint, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values = infect_clones, na.value=NA) +
  labs(x =  "cluster", y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Vaccine_imprints_VlnPlot_CloneType_infect.pdf",width = 5,height = 6)

z<-ggplot(subset@meta.data, aes(x=timepoint,y=Infection_imprint, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values = infect_clones, na.value=NA) +
  labs(x =  "cluster", y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z
ggsave(z,filename = "~/S3data/Images/Images/Infection_imprints_VlnPlot_CloneType_infect.pdf",width = 5,height = 6)

export <- subset@meta.data[,c("barcode","timepoint","Infection_imprint","highlight")]
export.wide <- export %>% unite("Clone_timepoint", highlight, timepoint, remove = T)
export.wide <- spread(export.wide, Clone_timepoint, Infection_imprint)
View(export.wide)
write.csv(export,"~/S3data/Export_csv/Infection_Top_Clones_GSVA_Scores.csv", row.names = F)


AIM.reactive.CD4@meta.data$highlight <- as.factor(ifelse(is.na(AIM.reactive.CD4@meta.data$highlight), "NA", AIM.reactive.CD4@meta.data$highlight))

z <- DimPlot(AIM.reactive.CD4, group.by = "highlight", split.by = "PriorCOVID",
             raster = FALSE,pt.size = 0.5,order = c("CAADFSGGGADGLTF_CASRTGSTNEKLFF", 
                                                                                    "CAMRAYNTGNQFYF_CASTSQGASDTEAFF",
                                                                                    "CAMRPLNTGNQFYF_CASSQESAGGIDEQFF",
                                                                                    "CAVYTSGTYKYIF_CSVDGQGNTGELFF", 
                                                                                    "CAVQDLLASGSRLTF_CASSSPSGELFF", 
                                                                                    "CATDWNNDMRF_CASSFTGTGNEKLFF",	
                                                                                    "CAMREVNTGTASKLTF_CASIPPGSAQGAFQPQHF",
                                                                                    "CAGRRGKLIF_CASSYEGPFGEQFF", 
                                                                                    "CAMRLGGSYIPTF_CASSPPRLQDNEQFF", 
                                                                                    "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF","NA")) + 
  theme(plot.title = element_blank())+  scale_color_manual(values = colorblind_vector(11)) 
z

View(AIM.reactive.CD4@meta.data$RNA_snn_res.0.4)


z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=PriorCOVID,y=Vaccine_imprint, color =RNA_snn_res.0.4 )) +
  theme_classic() +  geom_hline(yintercept = 0,size = 0.25)+geom_boxplot()+
  scale_color_brewer(palette = "Paired") + labs(x =  "cluster",
                                                     y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(timepoint~RNA_snn_res.0.4)
z

z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=PriorCOVID,y=Infection_imprint, color =RNA_snn_res.0.4 )) +
  theme_classic() +geom_hline(yintercept = 0,size = 0.25)+geom_boxplot()+
  scale_color_brewer(palette = "Paired") + labs(x =  "cluster",
                                                y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25)) + facet_grid(timepoint~RNA_snn_res.0.4)
z

#-----------GSVA ISG versus proliferative signature -- leading edge------------
#define proliferative and ISG signatures based off of leading edge from both pre and post GSEA comparisons 
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
prolif <- c(prolif_MS.pre$SYMBOL, prolif_MS.post$SYMBOL,prolif_G2M.pre$SYMBOL,prolif_G2M.post$SYMBOL)
prolif <- unique(prolif)

ISGa.pre <- subset(IFNa_pre,CORE.ENRICHMENT == "Yes")
ISGa.post <- subset(IFNa_post,CORE.ENRICHMENT == "Yes")
ISGg.pre <- subset(IFNg_pre,CORE.ENRICHMENT == "Yes")
ISGg.post <- subset(IFNg_post,CORE.ENRICHMENT == "Yes")
ISG <- c(ISGa.pre$SYMBOL, ISGa.post$SYMBOL,ISGg.pre$SYMBOL, ISGg.post$SYMBOL)
ISG <- unique(ISG)

AIM.reactive.CD4@meta.data <- AIM.reactive.CD4@meta.data %>% unite('COVID_timepoint',PriorCOVID, timepoint,  remove=FALSE)
DefaultAssay(AIM.reactive.CD4) <- "RNA"
#AIM.reactive.CD4.preBoost <- NormalizeData(AIM.reactive.CD4.preBoost)
ActivatedT.mxt.RNA <- as.data.frame(AIM.reactive.CD4@assays$RNA@scale.data) #36601 genes
#View(ActivatedT.mxt)
ActivatedT.mxt.RNA <- as.matrix(ActivatedT.mxt.RNA)
#View(ActivatedT.mxt.RNA)
hist(ActivatedT.mxt.RNA)

imprints <- list(prolif,ISG)
names(imprints) = c("Proliferative","ISG")
View(imprints)
gsva.es.all <- gsva(ActivatedT.mxt.RNA, imprints, verbose=TRUE, method="gsva")

saveRDS(gsva.es.all, file = "~/S3data/Saved_Objects/gsva_prolif_ISG_imprints.RDS") 
gsva.es.all <- readRDS("~/S3data/Saved_Objects/gsva_infection_vax_imprints.RDS")

dim(gsva.es.all)

metaData <- colnames(ActivatedT.mxt.RNA)
data_to_add <- rownames_to_column(AIM.reactive.CD4@meta.data ,"barcodes")
data_to_add <- data_to_add[,c("barcodes","subject","timepoint","PriorCOVID", "COVID_timepoint")]
colnames(data_to_add) <- c("barcode","subject","timepoint","PriorCOVID", "COVID_timepoint") 

metaData <- merge(metaData, data_to_add, by.x = 1, by.y = 1)
## Determine the number of cells per sample

colnames(metaData) <- c("barcode", "subject","timepoint","PriorCOVID", "COVID_timepoint") 
metaData <- metaData %>% arrange(subject)

gsva.es.all <- gsva.es.all[ , metaData$barcode]


metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("barcode")
cols <-  list(
  PriorCOVID = c(yes = "#b5b2f1", no = "#ffc471" ))

pathway_heatmap <- pheatmap::pheatmap(gsva.es.all,annotation_col = metaData, # Add metadata labels!
                                      #cutree_cols = 2,
                                      scale = "row",
                                      cluster_cols = T,
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6, # Shrink the pathway labels
                                      #annotation_colors	= cols, 
                                      color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

gsva.results <- t(gsva.es.all)
data <- merge(metaData, gsva.results, by.x = 0, by.y = 0) %>% as.data.frame()

ggplot(data, aes(x=subject,y=ISG, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Infection imprint") + facet_grid(rows="timepoint")

ggplot(data, aes(x=subject,y=Proliferative, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Vaccine imprint") + facet_grid(rows="timepoint")

data$timepoint <- factor(data$timepoint, levels = c("7moPostVax","1moPostBoost","PostBreakThru"))
data$COVID_timepoint <- factor(data$COVID_timepoint, levels = c("no_7moPostVax","no_1moPostBoost","no_PostBreakThru","yes_7moPostVax","yes_1moPostBoost"))
z <- ggplot(data, aes(x=ISG,y=Proliferative, color = COVID_timepoint)) +
  theme_classic() + geom_point(size=0.05) + geom_density2d(binwidth = 0.25,color="black", size = 0.25) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(pre_vax,post_vax,breakthru,pre_infect,post_infect)) + labs(x =  "Infection signature enrichment score",
                                                                                           y="Vaccine signature enrichment score") +  
  facet_grid(PriorCOVID ~ timepoint) + NoLegend() + theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z

z <- ggplot(data, aes(x=ISG,y=Proliferative, color = COVID_timepoint)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,  size = 0.15)+ geom_point(size=0.05) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(pre_vax,post_vax,breakthru,pre_infect,post_infect)) + labs(x =  "Infection signature enrichment score",
                                                                                           y="Vaccine signature enrichment score") +  
  facet_grid(subject ~ timepoint) + NoLegend()+ theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))
z


gsva.results.1 <- gsva.results %>% as.data.frame() %>% rownames_to_column( "barcodes")
AIM.reactive.CD4@meta.data  <- rownames_to_column(AIM.reactive.CD4@meta.data, "barcodes")
AIM.reactive.CD4@meta.data <- merge(AIM.reactive.CD4@meta.data, gsva.results.1, by = "barcodes")
AIM.reactive.CD4@meta.data <- column_to_rownames(AIM.reactive.CD4@meta.data, "barcodes")
View(AIM.reactive.CD4@meta.data)

AIM.reactive.CD4@meta.data$timepoint <- factor(AIM.reactive.CD4@meta.data$timepoint , levels = c("7moPostVax" , "1moPostBoost", "PostBreakThru"))
z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=ISG,y=Proliferative, color = RNA_snn_res.0.4)) +
  theme_classic()  +geom_point(size=0.25)+ geom_density2d(binwidth = 0.5,color="black", size = 0.25) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_brewer(palette = "Paired") + labs(x =  "Infection signature enrichment score",
                                                y="Vaccine signature enrichment score") +  
  facet_grid(RNA_snn_res.0.4 ~ COVID_timepoint) + theme(line=element_line(colour = "black", size = 0.25))

z


combined.AIM.CD4 <- expression2List(AIM.reactive.CD4, split.by = "subject_timepoint")
x<- compareClonotypes(combined.AIM.CD4,number = 1,  
                      cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
View(x) #look at top 20 clones 
top <- x #%>% head(20)
top_20_clones <- top$Clonotypes 
top_20_clones <- as.character(top_20_clones)
top_20_clones

seurat <- highlightClonotypes(AIM.reactive.CD4, 
                              cloneCall= "aa",
                              sequence = top_20_clones)
#View(seurat@meta.data)
z <- ggplot(seurat@meta.data, aes(x=ISG,y=Proliferative, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="grey", size = 0.25)+ geom_point(size=3) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = rev(colorblind_vector(123)), na.value = NA) + labs(x =  "Infection signature enrichment score",
                                                                                y="Vaccine signature enrichment score") +  
  facet_grid(PriorCOVID ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25)) + NoLegend()

z

CAMRLGGSYIPTF_CASSPPRLQDNEQFF
AIM.reactive.CD4 <- highlightClonotypes(AIM.reactive.CD4, 
                                        cloneCall= "aa", 
                                        sequence = c(#"CAADFSGGGADGLTF_CASRTGSTNEKLFF", 
                                          #"CAMRAYNTGNQFYF_CASTSQGASDTEAFF",
                                          #"CAMRPLNTGNQFYF_CASSQESAGGIDEQFF",
                                          #  "CAVYTSGTYKYIF_CSVDGQGNTGELFF", 
                                          "CAVQDLLASGSRLTF_CASSSPSGELFF","CALRDGGSNYKLTF_CAWSASPDTQYF", #Subject1
                                          "CALLANQAGTALIF_CASSESNRGAYGYTF","CAVQALLNAGNMLTF_CASGARTGEQFF", #Subject1
                                          "CAVVPLGNTGKLIF_CASSQDLTGANVLTF", #Subject4
                                          "CAGAEETSGSRLTF_CASSLSGGSNSPLHF", #Subject5
                                          "CAASEAAGNKLTF_CASSQVGGMGAKNIQYF", #Subject6
                                          #"CAMREVNTGTASKLTF_CASIPPGSAQGAFQPQHF",
                                          # "CAGRRGKLIF_CASSYEGPFGEQFF" ,
                                          "CAMRLGGSYIPTF_CASSPPRLQDNEQFF" , # Subject14
                                          "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF", # Subject14
                                          "CATDWNNDMRF_CASSFTGTGNEKLFF",   #Subject8 
                                          "CALNTGGFKTIF_CASSEAGGAGRVNEQFF"  #Subject11
                                        ))
AIM.reactive.CD4.vax <- subset(AIM.reactive.CD4, PriorCOVID == "no")
AIM.reactive.CD4.vax <- highlightClonotypes(AIM.reactive.CD4.vax, 
                                            cloneCall= "aa", 
                                            sequence = c("CALRDGGSNYKLTF_CAWSASPDTQYF",
                                                         "CAVQDLLASGSRLTF_CASSSPSGELFF","CALLANQAGTALIF_CASSESNRGAYGYTF",#"CAVQALLNAGNMLTF_CASGARTGEQFF","CVVHGNNFNKFYF_CSARAGEGNYGYTF",#Subject1
                                                         "CAVVPLGNTGKLIF_CASSQDLTGANVLTF" #Subject4 
                                                         # "CAASQGAQKLVF_CASSQARLAPQFF",      #Subject5
                                                         # "CAASEAAGNKLTF_CASSQVGGMGAKNIQYF", "CAAFISGNTPLVF_CSTGGRGGETQYF" #Subject7
                                            ))   
vax_clones = c( "#f6aa1c", "#f25c54", "#bf3100","#550527") #, "orange", "yellow2", "red", "red3","brown" )

z <- ggplot(AIM.reactive.CD4.vax@meta.data, aes(x=ISG,y=Proliferative, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="lightgrey", size = 0.25)+ geom_point(size=2) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values=vax_clones , na.value=NA) + labs(x =  "Infection signature enrichment score",
                                                             y="Vaccine signature enrichment score") +  
  facet_grid( ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25))+
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))
z

pre <- prop.table(table(AIM.reactive.CD4.vax$highlight, AIM.reactive.CD4.vax$RNA_snn_res.0.4), margin = 2)%>% as.data.frame()
pre$Freq <- pre$Freq*100
colnames(pre) <- c("Clone","Cluster", "%")
z <- ggplot(pre, aes(x=`Clone`, y=`%`, fill = Cluster, group = Cluster)) +
  geom_col(width = 0.75, colour = "white") +
  #geom_area(aes(x = c("no" = 1.25, "yes" = 1.75)[`PriorCOVID`]), 
  # position = "fill", colour = "white", alpha = 0.4,
  #outline.type = "both") + 
  theme_classic()  + scale_fill_brewer(palette = "Paired") +
  theme(panel.background=element_rect(fill="white"), 
        axis.text=element_text(colour="black",size=20), legend.text=element_text(colour="black",size=20),
        axis.title.y.left = element_text(size = 20))+labs(x="Prior COVID?", title = "Pre booster")+ RotatedAxis()
z



subset <- subset(AIM.reactive.CD4.vax, highlight == "CALRDGGSNYKLTF_CAWSASPDTQYF" & PriorCOVID == "no" |
                   highlight == "CAVQDLLASGSRLTF_CASSSPSGELFF"& PriorCOVID == "no" |highlight == "CALLANQAGTALIF_CASSESNRGAYGYTF" & PriorCOVID == "no"|
                   highlight == "CAVVPLGNTGKLIF_CASSQDLTGANVLTF" & PriorCOVID == "no" ) 

infection_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(ISG)) 

vaccine_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Proliferative)) 

ggplot(vaccine_score_vax,aes(x=timepoint,y=`sum(Proliferative)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of vaccine signature enrichment score")
ggplot(infection_score_vax,aes(x=timepoint,y=`sum(ISG)`, color = highlight, group = highlight)) + geom_point()+geom_line()+ theme_classic()+ 
  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of infection signature enrichment score")

z <- ggplot(subset@meta.data, aes(x=timepoint,y=Proliferative, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+ 
  labs(x =  "cluster", y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight")  + scale_y_continuous(limits = c(-0.6,0.6))
z

z<- ggplot(subset@meta.data, aes(x=timepoint,y=ISG, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+
  labs(x =  "cluster", y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z


export <- subset@meta.data[,c("barcode","timepoint","Infection_imprint","highlight")]
export.wide <- export %>% unite("Clone_timepoint", highlight, timepoint, remove = T)
export.wide <- spread(export.wide, Clone_timepoint, Infection_imprint)
View(export.wide)
write.csv(export,"~/S3data/Export_csv/Vaccine_Top_Clones_GSVA_Scores.csv", row.names = F)

#infection-primed 
AIM.reactive.CD4.infect<- subset(AIM.reactive.CD4, PriorCOVID == "yes")
AIM.reactive.CD4.infect <- highlightClonotypes(AIM.reactive.CD4.infect, 
                                               cloneCall= "aa", 
                                               sequence = c(#"CATDWNNDMRF_CASSFTGTGNEKLFF", #Subject8
                                                 "CAMRLGGSYIPTF_CASSPPRLQDNEQFF",#Subject14, Subject12
                                                 "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF", #"Subject14
                                                 "CALNTGGFKTIF_CASSEAGGAGRVNEQFF", 
                                                 "CAADFSGGGADGLTF_CASRTGSTNEKLFF" #Subject11
                                                 #"CAVRYSLNSGYSTLTF_CASSELLAGGPNEQYF" #Subject14, Subject8  
                                               ))
infect_clones = c( "#E53F71","#8338ec","#E090DF", "#541675")

z <- ggplot(AIM.reactive.CD4.infect@meta.data, aes(x=ISG,y=Proliferative, color = highlight)) +
  theme_classic() +  geom_density2d(binwidth = 0.25,color="lightgrey", size = 0.25)+ geom_point(size=0.5) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = infect_clones, na.value=NA) + labs(x =  "Infection signature enrichment score",
                                                                 y="Vaccine signature enrichment score") +  
  facet_grid( ~ timepoint) + theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z


subset <- subset(AIM.reactive.CD4.infect, #highlight ==   "CATDWNNDMRF_CASSFTGTGNEKLFF" & PriorCOVID == "yes"|   #Subject8
                 highlight ==    "CAMRLGGSYIPTF_CASSPPRLQDNEQFF" & PriorCOVID == "yes"|  # Subject14
                   highlight ==   "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF" & PriorCOVID == "yes" | # Subject14
                   highlight ==    "CALNTGGFKTIF_CASSEAGGAGRVNEQFF" & PriorCOVID == "yes"|
                   highlight ==    "CAADFSGGGADGLTF_CASRTGSTNEKLFF" & PriorCOVID == "yes"
                 # highlight ==    "CAVRYSLNSGYSTLTF_CASSELLAGGPNEQYF" & PriorCOVID == "yes" 
)

infection_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(ISG)) 

vaccine_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Proliferative)) 

ggplot(vaccine_score_vax,aes(x=timepoint,y=`sum(Proliferative)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values = infect_clones, na.value=NA) +labs(x="time point", y="Sum of Vaccine Imprint")
ggplot(infection_score_vax,aes(x=timepoint,y=`sum(ISG)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values = infect_clones, na.value=NA) +labs(x="time point", y="Sum of Infection Imprint")

z<-ggplot(subset@meta.data, aes(x=timepoint,y=Proliferative, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values = infect_clones, na.value=NA) +
  labs(x =  "cluster", y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z

z<-ggplot(subset@meta.data, aes(x=timepoint,y=ISG, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values = infect_clones, na.value=NA) +
  labs(x =  "cluster", y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z

AIM.reactive.CD4@meta.data$highlight <- as.factor(ifelse(is.na(AIM.reactive.CD4@meta.data$highlight), "NA", AIM.reactive.CD4@meta.data$highlight))

z <- DimPlot(AIM.reactive.CD4, group.by = "highlight", split.by = "PriorCOVID",
             raster = FALSE,pt.size = 0.5,order = c("CAADFSGGGADGLTF_CASRTGSTNEKLFF", 
                                                    "CAMRAYNTGNQFYF_CASTSQGASDTEAFF",
                                                    "CAMRPLNTGNQFYF_CASSQESAGGIDEQFF",
                                                    "CAVYTSGTYKYIF_CSVDGQGNTGELFF", 
                                                    "CAVQDLLASGSRLTF_CASSSPSGELFF", 
                                                    "CATDWNNDMRF_CASSFTGTGNEKLFF",	
                                                    "CAMREVNTGTASKLTF_CASIPPGSAQGAFQPQHF",
                                                    "CAGRRGKLIF_CASSYEGPFGEQFF", 
                                                    "CAMRLGGSYIPTF_CASSPPRLQDNEQFF", 
                                                    "CAQAEAQGGSEKLVF_CASSVGGGRNEKLFF","NA")) + 
  theme(plot.title = element_blank())+  scale_color_manual(values = colorblind_vector(11)) 
z

View(AIM.reactive.CD4@meta.data$RNA_snn_res.0.4)


z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=PriorCOVID,y=Proliferative_hallmark, color =RNA_snn_res.0.4 )) +
  theme_classic() +  geom_hline(yintercept = 0,size = 0.25)+geom_boxplot()+
  scale_color_brewer(palette = "Paired") + labs(x =  "cluster",
                                                y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(~RNA_snn_res.0.4)
z

z <- ggplot(AIM.reactive.CD4@meta.data, aes(x=PriorCOVID,y=Infection_imprint, color =RNA_snn_res.0.4 )) +
  theme_classic() +geom_hline(yintercept = 0,size = 0.25)+geom_boxplot()+
  scale_color_brewer(palette = "Paired") + labs(x =  "cluster",
                                                y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25)) + facet_grid(timepoint~RNA_snn_res.0.4)
z

#-----------GSVA ISG versus proliferative signature -- full hallmark gene sets------------------
#GSVA -- AIM Reactive CD4 enrichment scores for hallmark gene sets 
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
 
AIM.reactive.CD4@meta.data <- AIM.reactive.CD4@meta.data %>% unite('COVID_timepoint',PriorCOVID, timepoint,  remove=FALSE)
DefaultAssay(AIM.reactive.CD4) <- "RNA"
#AIM.reactive.CD4.preBoost <- NormalizeData(AIM.reactive.CD4.preBoost)
ActivatedT.mxt.RNA <- as.data.frame(AIM.reactive.CD4@assays$RNA@scale.data) #36601 genes
#View(ActivatedT.mxt)
ActivatedT.mxt.RNA <- as.matrix(ActivatedT.mxt.RNA)
#View(ActivatedT.mxt.RNA)
hist(ActivatedT.mxt.RNA)

imprints <- list(prolif,ISG)
names(imprints) = c("Proliferative","ISG")
View(imprints)
gsva.es.all <- gsva(ActivatedT.mxt.RNA, imprints, verbose=TRUE, method="gsva")

saveRDS(gsva.es.all, file = "~/S3data/Saved_Objects/gsva_prolif_ISG_imprints_hallmark_GeneSets.RDS") 
gsva.es.all <- readRDS("~/S3data/Saved_Objects/gsva_infection_vax_imprints.RDS")

dim(gsva.es.all)

metaData <- colnames(ActivatedT.mxt.RNA)
data_to_add <- rownames_to_column(AIM.reactive.CD4@meta.data ,"barcodes")
data_to_add <- data_to_add[,c("barcodes","subject","timepoint","PriorCOVID", "COVID_timepoint")]
colnames(data_to_add) <- c("barcode","subject","timepoint","PriorCOVID", "COVID_timepoint") 

metaData <- merge(metaData, data_to_add, by.x = 1, by.y = 1)
## Determine the number of cells per sample

colnames(metaData) <- c("barcode", "subject","timepoint","PriorCOVID", "COVID_timepoint") 
metaData <- metaData %>% arrange(subject)

gsva.es.all <- gsva.es.all[ , metaData$barcode]


metaData <- metaData %>%
  # pheatmap will want our sample names that match our data to
  tibble::column_to_rownames("barcode")
cols <-  list(
  PriorCOVID = c(yes = "#b5b2f1", no = "#ffc471" ))

pathway_heatmap <- pheatmap::pheatmap(gsva.es.all,annotation_col = metaData, # Add metadata labels!
                                      #cutree_cols = 2,
                                      scale = "row",
                                      cluster_cols = T,
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6, # Shrink the pathway labels
                                      #annotation_colors	= cols, 
                                      color = viridis(10))#colorRampPalette(c("blue", "white", "firebrick3"))(50))

gsva.results <- t(gsva.es.all)
data <- merge(metaData, gsva.results, by.x = 0, by.y = 0) %>% as.data.frame()

ggplot(data, aes(x=subject,y=ISG, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Infection imprint") + facet_grid(rows="timepoint")

ggplot(data, aes(x=subject,y=Proliferative, color = PriorCOVID)) + geom_hline(yintercept = 0) + geom_jitter() + geom_boxplot() +theme_classic() +
  scale_color_manual(values = c( "#ffc471","#b5b2f1")) + labs(x = "Subject",y="Vaccine imprint") + facet_grid(rows="timepoint")

data$timepoint <- factor(data$timepoint, levels = c("7moPostVax","1moPostBoost","PostBreakThru"))
data$COVID_timepoint <- factor(data$COVID_timepoint, levels = c("no_7moPostVax","no_1moPostBoost","no_PostBreakThru","yes_7moPostVax","yes_1moPostBoost"))
z <- ggplot(data, aes(x=ISG,y=Proliferative, color = COVID_timepoint)) +
  theme_classic() + geom_point(size=0.05) + geom_density2d(binwidth = 0.25,color="black", size = 0.25) +
  geom_hline(yintercept = 0,size = 0.25)+geom_vline(xintercept = 0, size = 0.25)+
  scale_color_manual(values = c(pre_vax,post_vax,breakthru,pre_infect,post_infect)) + labs(x =  "ISG enrichment score",
                                                                                           y="Proliferative enrichment score") +  
  facet_grid(PriorCOVID ~ timepoint) + NoLegend() + theme(line=element_line(colour = "black", size = 0.25)) +
  scale_y_continuous(limits = c(-0.6,0.6))+ scale_x_continuous(limits = c(-0.6,0.6))

z

colnames(gsva.results) <- c("Proliferative_hallmark","ISG_hallmark")
head(gsva.results)
gsva.results.1 <- gsva.results %>% as.data.frame() %>% rownames_to_column( "barcodes")
AIM.reactive.CD4@meta.data  <- rownames_to_column(AIM.reactive.CD4@meta.data, "barcodes")
AIM.reactive.CD4@meta.data <- merge(AIM.reactive.CD4@meta.data, gsva.results.1, by = "barcodes")
AIM.reactive.CD4@meta.data <- column_to_rownames(AIM.reactive.CD4@meta.data, "barcodes")
View(AIM.reactive.CD4@meta.data)

vax_clones = c( "#f6aa1c", "#f25c54", "#bf3100","#550527") #, "orange", "yellow2", "red", "red3","brown" )

subset <- subset(AIM.reactive.CD4, highlight == "CALRDGGSNYKLTF_CAWSASPDTQYF" & PriorCOVID == "no" |
                   highlight == "CAVQDLLASGSRLTF_CASSSPSGELFF"& PriorCOVID == "no" |highlight == "CALLANQAGTALIF_CASSESNRGAYGYTF" & PriorCOVID == "no"|
                   highlight == "CAVVPLGNTGKLIF_CASSQDLTGANVLTF" & PriorCOVID == "no" ) 

infection_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(ISG_hallmark)) 

vaccine_score_vax <- subset@meta.data %>% group_by(timepoint, highlight) %>%
  dplyr::summarize(sum(Proliferative_hallmark)) 

ggplot(vaccine_score_vax,aes(x=timepoint,y=`sum(Proliferative_hallmark)`, color = highlight, group = highlight)) + geom_point()+geom_line()+
  theme_classic()+  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of vaccine signature enrichment score")
ggplot(infection_score_vax,aes(x=timepoint,y=`sum(ISG)`, color = highlight, group = highlight)) + geom_point()+geom_line()+ theme_classic()+ 
  scale_color_manual(values=vax_clones , na.value=NA) +labs(x="time point", y="Sum of infection signature enrichment score")

z <- ggplot(subset@meta.data, aes(x=timepoint,y=Proliferative_hallmark, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+ 
  labs(x =  "cluster", y="Vaccine signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight")  + scale_y_continuous(limits = c(-0.6,0.6))
z

z<- ggplot(subset@meta.data, aes(x=timepoint,y=ISG_hallmark, color = highlight, group = highlight)) + geom_violin(aes(group = timepoint))+
  geom_jitter(width = 0.2, size = 0.75) +  theme_classic() +  
  geom_hline(yintercept = 0,size = 0.25,linetype = "dashed", color = "grey" )+
  scale_color_manual(values=vax_clones , na.value=NA)+
  labs(x =  "cluster", y="Infection signature enrichment score") +  
  theme(line=element_line(colour = "black", size = 0.25))+facet_grid(rows="highlight") + scale_y_continuous(limits = c(-0.6,0.6))
z


#----------------ScRpertoire------------------------------------------
DimPlot(AIM.integrated.T.clean.CD4, group.by = "cloneType") +
  scale_color_manual(values = rev(colorblind_vector(5)), na.value="grey") + 
  theme(plot.title = element_blank())
Idents(AIM.integrated.T.clean.CD4) <- "integrated_cluster.IDs"
occupiedscRepertoire(AIM.integrated.T.clean.CD4, x.axis = "integrated_cluster.IDs")

Idents(AIM.reactive.CD4) <- "subject_timepoint"
combined.AIM.CD4 <- expression2List(AIM.reactive.CD4, split.by = "subject_timepoint")
quantContig(combined.AIM.CD4, cloneCall="aa", chain = "TRB",
            scale = T)+NoLegend()
x <- clonalDiversity(combined.AIM.CD4, group.by = "subject",x.axis = "integrated_cluster.IDs",chain = "TRB",cloneCall = "aa",  exportTable = T, 
                n.boots = 500)
x
View(combined.AIM.CD4)
ggplot(x, aes(x = x$integrated_cluster.IDs, y = x$Shannon, fill = "subject")) + geom_boxplot() + theme_classic() + 
  scale_y_continuous(limits = c(1.5,5))

z<-clonalOverlap(combined.AIM.CD4, cloneCall="aa", chain = "TRB",method = "morisita")+
  labs(title="Clonal Overlap among AIM Reactive CD4+ T cluster")+RotatedAxis()
z
ggsave(z,filename = "./Images/Images//clonalOverlap_among_subjects_AIMreactiveCD4_cluster.pdf", width = 10, height = 10) 
combined.AIM.CD4

order <- c("Subject1_7moPostVax","Subject1_1moPostBoost","Subject1_PostBreakThru",
"Subject2_7moPostVax","Subject2_1moPostBoost" ,"Subject2_PostBreakThru", 
"Subject3_7moPostVax" ,"Subject3_1moPostBoost" ,
"Subject4_7moPostVax", "Subject4_1moPostBoost","Subject4_PostBreakThru",
"Subject5_7moPostVax" ,   "Subject5_1moPostBoost"  , "Subject5_PostBreakThru",
"Subject6_7moPostVax" , "Subject6_1moPostBoost" ,"Subject6_PostBreakThru", 
"Subject7_7moPostVax","Subject7_1moPostBoost" ,
"Subject8_7moPostVax"  ,  "Subject8_1moPostBoost",  
"Subject9_7moPostVax","Subject9_1moPostBoost",
"Subject10_7moPostVax"   ,"Subject10_1moPostBoost"  ,
"Subject11_7moPostVax"   ,"Subject11_1moPostBoost" ,
"Subject12_7moPostVax"   ,"Subject12_1moPostBoost" ,
"Subject13_7moPostVax"   ,   "Subject13_1moPostBoost" , 
"Subject14_7moPostVax" ,    "Subject14_1moPostBoost")
combined.AIM.CD4 <- combined.AIM.CD4[order]

Subject11 <- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject11_7moPostVax", "Subject11_1moPostBoost"), 
                  cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
Subject11$Sample <- factor(Subject11$Sample, levels = c("Subject11_7moPostVax", "Subject11_1moPostBoost"))
plot <- ggplot(Subject11, aes(x = Sample, fill = Clonotypes, 
                          group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                          y = Proportion, label = Clonotypes)) + theme_classic() + 
  theme(axis.title.x = element_blank())+ scale_fill_manual(values = plasma(15))
plot = plot + geom_stratum() + geom_flow(stat = "alluvium") + labs(title = "Subject 11")
plot


Subject14 <- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject14_7moPostVax", "Subject14_1moPostBoost"), 
                  cloneCall="aa", graph = "alluvial", chain = "TRB", exportTable = T)
Subject14$Sample <- factor(Subject14$Sample, levels = c("Subject14_7moPostVax","Subject14_1moPostBoost"))
plot.1 <- ggplot(Subject14, aes(x = Sample, fill = Clonotypes, 
                           group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                           y = Proportion, label = Clonotypes)) + theme_classic() + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = plasma(13))
plot.1 = plot.1 + geom_stratum() + geom_flow(stat = "alluvium") +labs(title = "Subject 14")
plot.1

Subject1 <- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject1_7moPostVax", "Subject1_1moPostBoost","Subject1_PostBreakThru" ), 
                  cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
Subject1$Sample <- factor(Subject1$Sample, levels = c("Subject1_7moPostVax", "Subject1_1moPostBoost","Subject1_PostBreakThru"))
plot.2 <- ggplot(Subject1, aes(x = Sample, fill = Clonotypes, 
                            group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                            y = Proportion, label = Clonotypes)) + theme_classic() + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = plasma(20))
plot.2 = plot.2 + geom_stratum() + geom_flow(stat = "alluvium") + labs(title = "Subject 1")
plot.2

Subject4<- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject4_7moPostVax", "Subject4_1moPostBoost","Subject4_PostBreakThru" ), 
                  cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
Subject4$Sample <- factor(Subject4$Sample, levels = c("Subject4_7moPostVax", "Subject4_1moPostBoost","Subject4_PostBreakThru"))
plot.3 <- ggplot(Subject4, aes(x = Sample, fill = Clonotypes, 
                           group = Clonotypes, stratum = Clonotypes, alluvium = Clonotypes, 
                           y = Proportion, label = Clonotypes)) + theme_classic() + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values = plasma(21))
plot.3 = plot.3 + geom_stratum() + geom_flow(stat = "alluvium") + labs(title = "Subject 4")
plot.3
View(Subject4)

ggsave(plot, filename = "./Images/Images/DominantClones_Subject11.pdf", width = 5, height = 5.25)
ggsave(plot.1, filename = "./Images/Images/DominantClones_Subject14.pdf", width = 5, height = 5.25)
ggsave(plot.2, filename = "./Images/Images/DominantClones_Subject1.pdf", width = 5, height = 5.25)
ggsave(plot.3, filename = "./Images/Images/DominantClones_Subject4.pdf", width = 5, height = 5.25)

Subject5 <- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject5_7moPostVax", "Subject5_1moPostBoost", "Subject5_PostBreakThru"), 
                           cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)
View(Subject5)
Subject6 <- compareClonotypes(combined.AIM.CD4, numbers = 10, samples = c("Subject6_7moPostVax", "Subject6_1moPostBoost", "Subject6_PostBreakThru"), 
                          cloneCall="aa", graph = "alluvial", chain = "both", exportTable = T)

z<- clonalProportion(combined.AIM.CD4, cloneCall = "aa",
                 split = c(1,5, 20, 100, 500), "TRB") + 
  theme(axis.text=element_text(colour="black",size=15))+ RotatedAxis()
z 
ggsave(z,filename = "./Images/Images/clonalProportion_among_subjects_AIMreactiveCD4_cluster.pdf", 
       width = 12, height = 8) 
