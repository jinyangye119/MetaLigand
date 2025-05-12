#' Score LR interactions for specified cl-pair in a specific sample condition.
#'
#' @param seuratobj Seurat object of the input scRNA-seq data, with assays "RNA" and "NPL_assay"; sample conditions and cell-types/clusters in the columns "cond" and "cl" of the meta.data.
#' @param lr_network a ligand_receptor network matrix with columns "L" and "R".
#' @param cl_from sender cell-type/cluster.
#' @param cl_to receiver cell-type/cluster.
#' @param cond a sample condition.
#' @param avgexprall the average expression of all genes in all cells.
#' @param LRI.method "cpdb", "natmi", or "scsigr".
#'
#' @return
#' LR interaction score data.frame for specified cell type/cluster-pair and sample condition
#' @import Seurat stringr tidyverse
#' @export

getLRIscorematrix <- function(seuratobj, lr_network, cl_from, cl_to, cond, avgexprall, LRI.method) {
  
  DefaultAssay(seuratobj) = "RNA"
  
  myobj = subset(seuratobj, subset = cl == cl_from)
  Idents(myobj) = 'cond'
  avg_expr_obj = AverageExpression(myobj, group.by = 'cond', return.seurat = TRUE)
  avg_expr_RNA = avg_expr_obj[['RNA']]$data[,cond,drop=FALSE]
  avg_expr_NPL = avg_expr_obj[['NPL_assay']]$data[,cond,drop=FALSE]
  avg_expr_from = rbind.data.frame(avg_expr_RNA, avg_expr_NPL)
  myobj = subset(seuratobj, subset = cl == cl_to)
  Idents(myobj) = 'cond'
  avg_expr_obj = AverageExpression(myobj, group.by = 'cond', return.seurat = TRUE)
  avg_expr_RNA = avg_expr_obj[['RNA']]$data[,cond,drop=FALSE]
  avg_expr_NPL = avg_expr_obj[['NPL_assay']]$data[,cond,drop=FALSE]
  avg_expr_to = rbind.data.frame(avg_expr_RNA, avg_expr_NPL)
  
  scorematrix = matrix(NA, nrow = nrow(lr_network), ncol = (ncol(lr_network) + 1))
  rownames(scorematrix) = sprintf("%s_%s", lr_network[,"L"], lr_network[,"R"])
  colnames(scorematrix) = c(colnames(lr_network), "SLRI")
  scorematrix = data.frame(scorematrix)
  for (i in 1:nrow(lr_network)) {
    scorematrix[i,colnames(lr_network)] = lr_network[i,]
    Lunits = str_split(lr_network[i,"L"], ";")[[1]]
    Rgenes = str_split(lr_network[i,"R"], ";")[[1]]
    if (length(intersect(Lunits, rownames(avg_expr_from))) == length(Lunits) &
        length(intersect(Rgenes, rownames(avg_expr_to))) == length(Rgenes)) {
      Lexpr = min(avg_expr_from[Lunits, cond])
      Rexpr = min(avg_expr_to[Rgenes, cond])
      if (LRI.method == "cpdb") {
        scorematrix[i,"SLRI"] = mean(c(Lexpr, Rexpr))
      } else if (LRI.method == "natmi") {
        scorematrix[i,"SLRI"] = Lexpr * Rexpr
      } else if (LRI.method == "scsigr") {
        scorematrix[i,"SLRI"] = sqrt(Lexpr * Rexpr)/(sqrt(Lexpr * Rexpr) + avgexprall)
      }
    }
  }
  scorematrix = scorematrix[!is.na(scorematrix$SLRI),]
  
  return(scorematrix)
}



