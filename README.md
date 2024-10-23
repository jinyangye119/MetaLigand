
# MetaLigand

MetaLigand is a full-featured R package for estimating non-peptide ligands abundance from bulk & single-cell RNA-seq data.
(MetaLigand is currently in beta and will be in active development through the peer review process.)

<hr>

<div  align="center">
<img src="Figures/Fig1.png" width = "850" height = "300" alt="LRLoops" align=center />
</div>

<hr>

## Quick Installation of LRLoop

```{r}
# Install devtools
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install required packages
install.packages(c('tidyverse','Seurat','ggplot2','pheatmap'))

# Install MetaLigand
devtools::install_github("https://github.com/jinyangye119/MetaLigand")
```

<hr>

## Learn to use MetaLigand

**load datasets and packages**
- List of NPLs and their associated synthetic and transporter genes were stored in .csv format under inst/extdata/
- List of datasets used in vignettes are availble under vignettes/

```{r}
library(tidyverse)
library(MetaLigand)
library(pheatmap)
library(Seurat)

```

**Testing MetaLigand on a single-cell RNA-seq dataset**
- A single-cell transcriptomic database of adult mouse visual cortex was used for testing (https://doi.org/10.1038/s41586-018-0654-5)
- Raw gene expression matrix and metadata were downloaded from GEO-GSE115746 and stored under vignettes/

```{r}
# Load VISP dataset
visp <- read.csv("vignettes/GSE115746_cells_exon_counts.csv.gz",header = T,row.names = 1)
visp_meta <- read.csv("vignettes/GSE115746_complete_metadata_28706-cells.csv.gz")%>%
  dplyr::filter(sample_name%in%colnames(Visp))%>%
  column_to_rownames("sample_name")
Visp <- Visp[,colnames(Visp)%in%rownames(meta_visp)]
```
```{r}
# This is standard Seurat pipeline, to visulize data before MetaLigand
# If you just want a NPL table, skip this step
Seurat<- CreateSeuratObject(
  counts = Visp,
  meta.data = meta_visp
)
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat)
Seurat <- ScaleData(Seurat)
Seurat <- RunPCA(Seurat)
Seurat <- FindNeighbors(Seurat,reduction = "pca", dims = 1:30)
Seurat <- FindClusters(Seurat,reduction = "pca", resolution = 0.5)
Seurat <- RunUMAP(Seurat,reduction = "pca", dims = 1:30)
Seurat <- subset(Seurat,cell_subclass%in%c("","Batch Grouping","Doublet","High Intronic","No Class","Low Quality"),invert=T)
Seurat$new_group <- paste(Seurat$cell_class,Seurat$cell_subclass,sep = "_")
DimPlot(Seurat,group.by = "new_group",label = T)
```


# Expression based on single cell
ave_expr = Seurat@assays$RNA@data%>%
  as.matrix()

Singlecell_matrix = Mega_matrix(ave_expr,species = "mouse")
Singlecell_matrix = as.data.frame(Singlecell_matrix)%>%
  drop_na()
Singlecell_matrix = Singlecell_matrix[rowSums(Singlecell_matrix)>0,]%>%
  as.matrix()

Seurat_new = Seurat
NPL_assay = CreateAssayObject(data=Singlecell_matrix)

Seurat_new[["NPL_assay"]] = NPL_assay

Seurat_new$new_group = factor(Seurat_new$new_group)
Seurat_new$new_group_num = as.numeric(Seurat_new$new_group)
Seurat_new$new_group_anno = paste(Seurat_new$new_group_num,Seurat_new$new_group,sep = ": ")
Seurat_new$new_group_anno = factor(Seurat_new$new_group_anno,levels = paste(1:23,levels(Seurat_new$new_group),sep = ": "))
Seurat_new$new_group_num_factor = factor(Seurat_new$new_group_num,levels=c(1:23))

Idents(Seurat_new)<-"new_group_num"
DefaultAssay(Seurat_new) <- "RNA"
DimPlot(Seurat_new,group.by = "new_group_num",label = T)&
  ggplot2::theme(legend.position = "right",
                 panel.background  = element_rect(colour = "black", size=0.8),
                 text=element_text(size = 14),
                 axis.text = element_text(size = 14),
                 strip.text =element_text(size = 14),
                 legend.text = element_text(size=13,face = "plain"))&
  guides(color = guide_legend(ncol=2,override.aes=list(shape = 19,size=4)))&
  scale_color_hue(label=levels(Seurat_new$new_group_anno))&
  ggtitle("")

DefaultAssay(Seurat_new) <- "NPL_assay"
FeaturePlot(Seurat_new,features = c("L-Glutamic acid","gamma-Aminobutyric acid"),label = F,min.cutoff = "q10",ncol = 1)&
  ggplot2::theme(legend.position = "right",
                 panel.background  = element_rect(colour = "black", size=0.8),
                 text=element_text(size = 14),
                 axis.text = element_text(size = 14),
                 strip.text =element_text(size = 14),
                 legend.text = element_text(size=14,face = "plain"))

VlnPlot(Seurat_new,group.by = "new_group_num_factor",features = c("L-Glutamic acid","gamma-Aminobutyric acid"),pt.size=0,ncol = 1)&
  ggplot2::theme(legend.position = "right",
                 panel.background  = element_rect(colour = "black", size=0.8),
                 text=element_text(size = 14),
                 axis.text.x = element_text(size = 14,angle = 0,hjust = 0.5,vjust = 1),
                 axis.text.y = element_text(size = 14,angle = 0),
                 strip.text =element_text(size = 14),
                 legend.text = element_text(size=14,face = "plain"))&
  NoLegend()&
  xlab("")&
  ggtitle("")

# Re-clustering 
Seurat_new <- FindVariableFeatures(Seurat_new,assay = "NPL_assay")
Seurat_new <- ScaleData(Seurat_new,assay = "NPL_assay")
Seurat_new <- RunPCA(Seurat_new,assay = "NPL_assay")
Seurat_new <- FindNeighbors(Seurat_new,assay = "NPL_assay",reduction = "pca", dims = 1:10)
Seurat_new <- FindClusters(Seurat_new,assay="NPL_assay",reduction = "pca", resolution = 0.1)
Seurat_new <- RunUMAP(Seurat_new,assay="NPL_assay",reduction = "pca", dims = 1:30)

DefaultAssay(Seurat_new) <- "NPL_assay"
DimPlot(Seurat_new,group.by = "new_group_num",label = T)&
  ggplot2::theme(legend.position = "right",
                 panel.background  = element_rect(colour = "black", size=0.8),
                 text=element_text(size = 14),
                 axis.text = element_text(size = 14),
                 strip.text =element_text(size = 14),
                 legend.text = element_text(size=13,face = "plain"))&
  guides(color = guide_legend(ncol=2,override.aes=list(shape = 19,size=4)))&
  scale_color_hue(label=levels(Seurat_new$new_group_anno))&
  ggtitle("")

FeaturePlot(Seurat_new,features = c("L-Glutamic acid","gamma-Aminobutyric acid"),label = F,min.cutoff = "q10",ncol = 1)&
  ggplot2::theme(legend.position = "right",
                 panel.background  = element_rect(colour = "black", size=0.8),
                 text=element_text(size = 14),
                 axis.text = element_text(size = 14),
                 strip.text =element_text(size = 14),
                 legend.text = element_text(size=14,face = "plain"))


```




<hr>

## Issues using LRLoop?

LRLoop is currently in __beta__. If you think you have found a bug, please [report an issue on Github](https://github.com/Pinlyu3/LRLoop/issues) with the __Bug Report__ form.
