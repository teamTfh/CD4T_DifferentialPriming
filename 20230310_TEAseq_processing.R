library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)
#BiocManager::install("biovizBase")
library(biovizBase)
library(EnsDb.Hsapiens.v86)
#BiocManager::install("Bioconductor/GenomeInfoDb", force = T)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
#BiocManager::install(version = "3.16")
library(BiocManager)
#install.packages("data.table")
library(data.table)
#BiocManager::install("scDblFinder")
library(scDblFinder)
library(tibble)
BiocManager::install("glmGamPoi")
library(glmGamPoi)

sessionInfo()

setwd("~/S3data/TEAseq")
#----------creating objects---------------
stim1.out <- Read10X(data.dir = "./20230207_Stim1/filtered_feature_bc_matrix")
stim2.out <- Read10X(data.dir = "./20230207_Stim2/filtered_feature_bc_matrix")
stim3.out <- Read10X(data.dir = "./20230207_Stim3/filtered_feature_bc_matrix")
unstim.out <- Read10X(data.dir = "./20230207_Unstim/filtered_feature_bc_matrix")

stim1 <- CreateSeuratObject(stim1.out$`Gene Expression`)
stim2 <- CreateSeuratObject(stim2.out$`Gene Expression`)
stim3 <- CreateSeuratObject(stim3.out$`Gene Expression`)
unstim <- CreateSeuratObject(unstim.out$`Gene Expression`)

#----------adding ADT assay-----------
Stim1_adt_mat <- readRDS("./20230207_Stim1/Stim1_adt_count_mat.rds")
Stim1_adt_mat <- Matrix::t(as(Stim1_adt_mat, "dgCMatrix"))
Stim1_adt_mat <- Stim1_adt_mat[!grepl("Control", rownames(Stim1_adt_mat)),]
stim1[["ADT"]] <- CreateAssayObject(counts = Stim1_adt_mat)

Stim2_adt_mat <- readRDS("./20230207_Stim2/Stim2_adt_count_mat.rds")
Stim2_adt_mat <- Matrix::t(as(Stim2_adt_mat, "dgCMatrix"))
Stim2_adt_mat <- Stim2_adt_mat[!grepl("Control", rownames(Stim2_adt_mat)),]
stim2[["ADT"]] <- CreateAssayObject(counts = Stim2_adt_mat)

Stim3_adt_mat <- readRDS("./20230207_Stim3/Stim3_adt_count_mat.rds")
Stim3_adt_mat <- Matrix::t(as(Stim3_adt_mat, "dgCMatrix"))
Stim3_adt_mat <- Stim3_adt_mat[!grepl("Control", rownames(Stim3_adt_mat)),]
stim3[["ADT"]] <- CreateAssayObject(counts = Stim3_adt_mat)

Unstim_adt_mat <- readRDS("./20230207_Unstim/Unstim_adt_count_mat.rds")
Unstim_adt_mat <- Matrix::t(as(Unstim_adt_mat, "dgCMatrix"))
Unstim_adt_mat <- Unstim_adt_mat[!grepl("Control", rownames(Unstim_adt_mat)),]
unstim[["ADT"]] <- CreateAssayObject(counts = Unstim_adt_mat)

#----------adding HTO assay-----------
Stim1_HTO_mat <- readRDS("./20230207_Stim1/Stim1_HTO_count_mat.rds")
Stim1_HTO_mat <- Matrix::t(as(Stim1_HTO_mat, "dgCMatrix"))
Stim1_HTO_mat <- Stim1_HTO_mat[!grepl("Control", rownames(Stim1_HTO_mat)),]
stim1[["HTO"]] <- CreateAssayObject(counts = Stim1_HTO_mat)

Stim2_HTO_mat <- readRDS("./20230207_Stim2/Stim2_HTO_count_mat.rds")
Stim2_HTO_mat <- Matrix::t(as(Stim2_HTO_mat, "dgCMatrix"))
Stim2_HTO_mat <- Stim2_HTO_mat[!grepl("Control", rownames(Stim2_HTO_mat)),]
stim2[["HTO"]] <- CreateAssayObject(counts = Stim2_HTO_mat)

Stim3_HTO_mat <- readRDS("./20230207_Stim3/Stim3_HTO_count_mat.rds")
Stim3_HTO_mat <- Matrix::t(as(Stim3_HTO_mat, "dgCMatrix"))
Stim3_HTO_mat <- Stim3_HTO_mat[!grepl("Control", rownames(Stim3_HTO_mat)),]
stim3[["HTO"]] <- CreateAssayObject(counts = Stim3_HTO_mat)

Unstim_HTO_mat <- readRDS("./20230207_Unstim/Unstim_HTO_count_mat.rds")
Unstim_HTO_mat <- Matrix::t(as(Unstim_HTO_mat, "dgCMatrix"))
Unstim_HTO_mat <- Unstim_HTO_mat[!grepl("Control", rownames(Unstim_HTO_mat)),]
#remove HTOs not used in this lane 
Unstim_HTO_mat <- Unstim_HTO_mat[setdiff(rownames(x = Unstim_HTO_mat),c( "Subject1_PostSecondDose","Subject17_PostSecondDose", 
                                                                         ".Subject7_PostSecondDose",
                                                                         "Subject23_PostInfection","Subject10_PostInfection" )), ]
head(Unstim_HTO_mat)
unstim[["HTO"]] <- CreateAssayObject(counts = Unstim_HTO_mat)

#----------adding ATAC assay-----------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

fragpath.stim1 <-"./20230207_Stim1/atac_fragments.tsv.gz"
stim1[["ATAC"]] <- CreateChromatinAssay(
  counts = stim1.out$Peaks,
  sep = c(":", "-"),
  fragments = fragpath.stim1,
  annotation = annotation
)

fragpath.stim2 <-"./20230207_Stim2/atac_fragments.tsv.gz"
stim2[["ATAC"]] <- CreateChromatinAssay(
  counts = stim2.out$Peaks,
  sep = c(":", "-"),
  fragments = fragpath.stim2,
  annotation = annotation
)

fragpath.stim3 <-"./20230207_Stim3/atac_fragments.tsv.gz"
stim3[["ATAC"]] <- CreateChromatinAssay(
  counts = stim3.out$Peaks,
  sep = c(":", "-"),
  fragments = fragpath.stim3,
  annotation = annotation
)

fragpath.unstim <-"./20230207_Unstim/atac_fragments.tsv.gz"
unstim[["ATAC"]] <- CreateChromatinAssay(
  counts = unstim.out$Peaks,
  sep = c(":", "-"),
  fragments = fragpath.unstim,
  annotation = annotation
)

#-------HTODemux--------
DefaultAssay(stim1) <- "HTO"
stim1 <- NormalizeData(stim1, normalization.method = "CLR", margin = 2)
stim1 <- HTODemux(stim1, assay = "HTO", positive.quantile = 0.99)

DefaultAssay(stim2) <- "HTO"
stim2 <- NormalizeData(stim2, normalization.method = "CLR", margin = 2)
stim2 <- HTODemux(stim2, assay = "HTO", positive.quantile = 0.99)

DefaultAssay(stim3) <- "HTO"
stim3 <- NormalizeData(stim3, normalization.method = "CLR", margin = 2)
stim3 <- HTODemux(stim3, assay = "HTO", positive.quantile = 0.99)

DefaultAssay(unstim) <- "HTO"
unstim <- NormalizeData(unstim, normalization.method = "CLR", margin = 2)
rowSums(unstim[["HTO"]])
unstim <- HTODemux(unstim, assay = "HTO", positive.quantile = 0.99)

## Determine how less stringent demultiplex would look with HTO_maxID at or above a set UMI threshold of ??????
stim1@meta.data <- stim1@meta.data  %>%
  mutate(HTO_classification.global = case_when(nCount_HTO >= 10 & HTO_classification.global == "Negative"  ~ "Singlet",
    nCount_HTO < 10 & HTO_classification.global == "Negative"  ~ "Negative",HTO_classification.global == "Singlet"  ~ "Singlet",
    HTO_classification.global == "Doublet"  ~ "Doublet"))

stim2@meta.data <- stim2@meta.data  %>%
  mutate(HTO_classification.global = case_when(nCount_HTO >= 10 & HTO_classification.global == "Negative"  ~ "Singlet",
                                               nCount_HTO < 10 & HTO_classification.global == "Negative"  ~ "Negative",HTO_classification.global == "Singlet"  ~ "Singlet",
                                               HTO_classification.global == "Doublet"  ~ "Doublet"))

stim3@meta.data <- stim3@meta.data  %>%
  mutate(HTO_classification.global = case_when(nCount_HTO >= 10 & HTO_classification.global == "Negative"  ~ "Singlet",
                                               nCount_HTO < 10 & HTO_classification.global == "Negative"  ~ "Negative",HTO_classification.global == "Singlet"  ~ "Singlet",
                                               HTO_classification.global == "Doublet"  ~ "Doublet"))

unstim@meta.data <- unstim@meta.data  %>%
  mutate(HTO_classification.global = case_when(nCount_HTO >= 10 & HTO_classification.global == "Negative"  ~ "Singlet",
                                               nCount_HTO < 10 & HTO_classification.global == "Negative"  ~ "Negative",HTO_classification.global == "Singlet"  ~ "Singlet",
                                               HTO_classification.global == "Doublet"  ~ "Doublet"))

#-----------------------------------------filter out doublets--------------------------------------------------
# using scDblFinder #doublet rate for multiome kit...assumed to be 0.8% per 1,000 cells recovered 
DefaultAssay(stim1) <- "RNA"
sce.combined <- as.SingleCellExperiment(stim1)
sce.combined <- scDblFinder(sce.combined, dbr = (0.008 * ncol(sce.combined)/1000), samples="HTO_maxID")
table(sce.combined$scDblFinder.class, sce.combined$HTO_maxID)
identical(colnames(stim1),colnames(sce.combined))
stim1$scDblFinder.class <- sce.combined$scDblFinder.class
stim1$scDblFinder.class <- factor(stim1$scDblFinder.class, levels=c("doublet", "singlet"))

DefaultAssay(stim2) <- "RNA"
sce.combined <- as.SingleCellExperiment(stim2)
sce.combined <- scDblFinder(sce.combined, dbr = (0.008 * ncol(sce.combined)/1000), samples="HTO_maxID")
table(sce.combined$scDblFinder.class, sce.combined$HTO_maxID)
identical(colnames(stim2),colnames(sce.combined))
stim2$scDblFinder.class <- sce.combined$scDblFinder.class
stim2$scDblFinder.class <- factor(stim2$scDblFinder.class, levels=c("doublet", "singlet"))

DefaultAssay(stim3) <- "RNA"
sce.combined <- as.SingleCellExperiment(stim3)
sce.combined <- scDblFinder(sce.combined, dbr = (0.008 * ncol(sce.combined)/1000), samples="HTO_maxID")
table(sce.combined$scDblFinder.class, sce.combined$HTO_maxID)
identical(colnames(stim3),colnames(sce.combined))
stim3$scDblFinder.class <- sce.combined$scDblFinder.class
stim3$scDblFinder.class <- factor(stim3$scDblFinder.class, levels=c("doublet", "singlet"))

DefaultAssay(unstim) <- "RNA"
sce.combined <- as.SingleCellExperiment(unstim)
sce.combined <- scDblFinder(sce.combined, dbr = (0.008 * ncol(sce.combined)/1000), samples="HTO_maxID")
table(sce.combined$scDblFinder.class, sce.combined$HTO_maxID)
identical(colnames(unstim),colnames(sce.combined))
unstim$scDblFinder.class <- sce.combined$scDblFinder.class
unstim$scDblFinder.class <- factor(unstim$scDblFinder.class, levels=c("doublet", "singlet"))

#----rename and save objects--------
obj_stim1_2 <- stim1
obj_stim2_2 <- stim2
obj_stim3_2 <- stim3
obj_unstim_2 <- unstim

saveRDS(obj_stim1_2, file = "obj_stim1_2.rds")
saveRDS(obj_stim2_2, file = "obj_stim2_2.rds")
saveRDS(obj_stim3_2, file = "obj_stim3_2.rds")
saveRDS(obj_unstim_2, file = "obj_unstim_2.rds")



