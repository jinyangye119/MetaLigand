#' Score LR interactions 
#'
#' @param lr_network a ligand_receptor network matrix with columns "L" and "R".
#' @param exprinfo output of function "getExpInfo".
#' @param LRI.method "cpdb", "natmi", or "scsigr".
#'
#' @return
#' List of LR interaction score data.frames for each cell type/cluster-pair in each sample condition
#' @import Seurat stringr tidyverse
#' @export


getLRIScore <- function(lr_network, exprinfo, LRI.method = "scsigr") {
  conds = colnames(exprinfo$avgexprlist[[1]])
  cls = names(exprinfo$avgexprlist)
  LRIscorelist = list()
  clpairidx = 0
  clpairids = vector()
  for (m in 1:length(cls)) {
    for(n in 1:length(cls)) {
      clpairidx = clpairidx + 1
      cl_from = cls[m]
      cl_to = cls[n]
      clpairids[clpairidx] = sprintf("%s_to_%s", cl_from, cl_to)
      print(sprintf("%s of %s, %s_to_%s", clpairidx, length(exprinfo$avgexprlist)^2, cl_from, cl_to))
      scorematrix = matrix(NA, nrow = nrow(lr_network), ncol = (ncol(lr_network) + length(conds)))
      rownames(scorematrix) = sprintf("%s_%s", lr_network[,"L"], lr_network[,"R"])
      colnames(scorematrix) = c(colnames(lr_network), sprintf("SLRI_%s", conds))
      scorematrix = data.frame(scorematrix)
      for (i in 1:nrow(lr_network)) {
        scorematrix[i,colnames(lr_network)] = lr_network[i,]
        Lunits = str_split(lr_network[i,"L"], ";")[[1]]
        Rgenes = str_split(lr_network[i,"R"], ";")[[1]]
        if (length(intersect(Lunits, rownames(exprinfo$avgexprlist[[1]]))) == length(Lunits) &
            length(intersect(Rgenes, rownames(exprinfo$avgexprlist[[1]]))) == length(Rgenes)) {
          for (j in 1:length(conds)) {
            Lexpr = min(exprinfo$avgexprlist[[cl_from]][Lunits, conds[j]])
            Rexpr = min(exprinfo$avgexprlist[[cl_to]][Rgenes, conds[j]])
            if (LRI.method == "cpdb") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = mean(c(Lexpr, Rexpr))
            } else if (LRI.method == "natmi") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = Lexpr * Rexpr
            } else if (LRI.method == "scsigr") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = sqrt(Lexpr * Rexpr)/(sqrt(Lexpr * Rexpr) + exprinfo$avgexprall)
            }
          }
        }
      }
      LRIscorelist[[clpairidx]] = scorematrix
    }
  }
  names(LRIscorelist) = clpairids
  
  return(LRIscorelist)
}




