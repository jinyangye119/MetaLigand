#' Get detection rates and average gene and NPL levels.
#'
#' @param seuratobj Seurat object of the input scRNA-seq data, with assays "RNA" and "NPL_assay"; sample conditions and cell-types/clusters in the columns "cond" and "cl" of the meta.data.
#'
#' @return
#' A list exprinfo:
#' exprinfo$avgexprall: the average expression of all genes in all cells
#' exprinfo$avgexprlist: list of average gene and NPL levels matrix of each cell-type/cluser in cls
#' exprinfo$pctexprlist: list of gene and NPL detection rate matrix of each cell-type/cluser in cls
#' @import Seurat tidyverse
#' @export


getExprInfo <- function(seuratobj) {
  cls = unique(seuratobj@meta.data$cl)
  conds = unique(seuratobj@meta.data$cond)
  avgall = log1p(mean(expm1(as.matrix(seuratobj[['RNA']]$data))))
  DefaultAssay(seuratobj) = "RNA"
  avgexprlist = list()
  pctexprlist = list()
  for (k in 1:length(cls)) {
    print(sprintf("k = %s (in 1:%s), cl = %s", k, length(cls), cls[k]))
    myobj = subset(seuratobj, subset = cl == cls[k])
    Idents(myobj) = 'cond'
    avg_expr_obj = AverageExpression(myobj, group.by = 'cond', return.seurat = TRUE)
    avg_expr_RNA = avg_expr_obj[['RNA']]$data[,conds,drop=FALSE]
    avg_expr_NPL = avg_expr_obj[['NPL_assay']]$data[,conds,drop=FALSE]
    avg_expr = rbind.data.frame(avg_expr_RNA, avg_expr_NPL)
    normdata = myobj[['RNA']]$data
    pct_expr = matrix(0, nrow = nrow(normdata), ncol = length(conds))
    rownames(pct_expr) = rownames(normdata)
    colnames(pct_expr) = conds
    for (i in 1:length(conds)) {
      subnormdata = normdata[,WhichCells(myobj, idents = conds[i]), drop=F]
      if (ncol(subnormdata) > 0) {
        pct_expr[,i] = rowSums(subnormdata>0)/ncol(subnormdata)
      }
    }
    pct_expr_RNA = pct_expr
    normdata = myobj[['NPL_assay']]$data
    pct_expr = matrix(0, nrow = nrow(normdata), ncol = length(conds))
    rownames(pct_expr) = rownames(normdata)
    colnames(pct_expr) = conds
    for (i in 1:length(conds)) {
      subnormdata = normdata[,WhichCells(myobj, idents = conds[i]), drop=F]
      if (ncol(subnormdata) > 0) {
        pct_expr[,i] = rowSums(subnormdata>0)/ncol(subnormdata)
      } 
    }
    pct_expr_NPL = pct_expr
    pct_expr = rbind.data.frame(pct_expr_RNA, pct_expr_NPL)
    avgexprlist[[k]] = avg_expr
    pctexprlist[[k]] = pct_expr
  }
  names(avgexprlist) = cls
  names(pctexprlist) = cls
  exprinfo = list(avgall, avgexprlist, pctexprlist)
  names(exprinfo) = c("avgexprall", "avgexprlist", "pctexprlist")

  return(exprinfo)
}




