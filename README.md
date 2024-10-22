
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
- List of datasets used in vignettes are availble under 

```{r}
library(tidyverse)
library(MetaLigand)
library(pheatmap)
library(Seurat)






```

**Download datasets used in vignettes**

To learn how to use LRLoop, read the following vignettes explaining several types of analyses:
- The datasets in vignettes can be downloaded from [google drive](https://drive.google.com/drive/folders/1WV0iSlAXCUwSZMSBnzsHdZc26RuyfunC?usp=sharing)

- [Perform LRLoop analysis starting from Seurat objects](vignettes/Main.md)





<hr>

## Issues using LRLoop?

LRLoop is currently in __beta__. If you think you have found a bug, please [report an issue on Github](https://github.com/Pinlyu3/LRLoop/issues) with the __Bug Report__ form.
