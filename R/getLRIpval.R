#' Score LR interactions for specified cl-pair in a specific sample condition, together with p-values.
#'
#' @param seuratobj Seurat object of the input scRNA-seq data, with assays "RNA" and "NPL_assay"; sample conditions and cell-types/clusters in the columns "cond" and "cl" of the meta.data.
#' @param lr_network a ligand_receptor network matrix with columns "L" and "R".
#' @param cl_from sender cell-type/cluster.
#' @param cl_to receiver cell-type/cluster.
#' @param cond a sample condition.
#' @param avgexprall the average expression of all genes in all cells.
#' @param LRI.method "cpdb", "natmi", or "scsigr".
#' @param numperm number of permutations.
#'
#' @return
#' LR interaction score and p-value data.frame for specified cell type/cluster-pair and sample condition
#' @import Seurat stringr tidyverse
#' @export


getLRIpval <- function(seuratobj, lr_network, cl_from, cl_to, cond, avgexprall, LRI.method = "scsigr", numperm = 100) {
  
  myscorematrix = getLRIscorematrix(seuratobj = seuratobj, lr_network = lr_network, 
                                    cl_from = cl_from, cl_to = cl_to, cond = cond, 
                                    avgexprall = avgexprall, LRI.method = LRI.method)
  
  condcells = rownames(seuratobj@meta.data)[seuratobj@meta.data$cond == cond]
  clvec = seuratobj@meta.data[condcells,"cl"]
  clvec = as.vector(clvec)
  names(clvec) = condcells
  
  permcls = vector()
  for (k in 1:numperm) {
    clvec_perm = sample(clvec, replace = F)
    names(clvec_perm) = names(clvec)
    permcls = cbind(permcls, clvec_perm)
  }
  rownames(permcls) = condcells
  colnames(permcls) = NULL
  
  permscores = matrix(0, nrow = nrow(myscorematrix), ncol = numperm)
  rownames(permscores) = rownames(myscorematrix)
  permscores = data.frame(permscores)
  for (k in 1:numperm) {
    print(sprintf("permutation No. %s", k))
    clvec_perm = permcls[,k]
    names(clvec_perm) = condcells
    seuratobjperm = seuratobj
    seuratobjperm@meta.data[names(clvec_perm),"cl"] = clvec_perm
    myscorematrixperm = getLRIscorematrix(seuratobj = seuratobjperm, lr_network = lr_network, 
                                          cl_from = cl_from, cl_to = cl_to, cond = cond, 
                                          avgexprall = avgexprall, LRI.method = LRI.method)
    permscores[,k] = myscorematrixperm[rownames(myscorematrix), "SLRI"]
  }
  
  myscorematrix$pval = NA
  for (i in 1:nrow(myscorematrix)) {
    myscorematrix[i,"pval"] = sum(permscores[i,] > myscorematrix[i,"SLRI"])/numperm
  }
  
  return(myscorematrix)
}


