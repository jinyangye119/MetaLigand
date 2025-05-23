#' Creat a ligand expression matrix
#'
#' This function loads gene expression matrix and calulate ligand expression
#' matrix based on gene expression. It assume that each columns are variables
#' and column names are variable names. It assumes that rownames are gene symbols.
#'
#' @param infile "./inst/extdata"
#' @return A matrix of the infile
#' @export

Meta_matrix <- function(ave_expr,
                        Synthetic_step2 = T,
                        precusor_transporter = T,
                        species = NULL){
  suppressMessages(require(tidyverse))

  `%notin%` <- Negate(`%in%`)

  gm_mean <- function(x){
    prod(x)^(1/length(x))
  }

  `%notin%` <- Negate(`%in%`)

  if (species%notin%c("mouse","human","zebrafish")){
    stop("species is required: mouse, human, zebrafish")
  }

  # Input files
  if (species == "mouse"){
    LR_data <- read.csv(system.file("extdata", "LR_data_mouse.csv", package = "MetaLigand"))
  } else if (species =="human"){
    LR_data <- read.csv(system.file("extdata", "LR_data_human.csv", package = "MetaLigand"))
  } else {
    LR_data <- read.csv(system.file("extdata", "LR_data_zebrafish.csv", package = "MetaLigand"))
  }

  # Creat all databases
  All_synthetic = LR_data%>%
    dplyr::filter(!Synthetic_genes=="")%>%
    dplyr::select(Compounds, Synthetic_genes)

  All_synthetic = split(All_synthetic,list(All_synthetic$Compounds))

  All_synthetic_lv2 = LR_data%>%
    dplyr::filter(!Synthetic_genes_lv2=="")%>%
    dplyr::select(Compounds, Synthetic_genes_lv2)

  All_synthetic_lv2 = split(All_synthetic_lv2,list(All_synthetic_lv2$Compounds))

  All_trans = LR_data%>%
    dplyr::filter(!Transporter_genes=="")%>%
    dplyr::select(Compounds, Transporter_genes)
  All_trans = split(All_trans,list(All_trans$Compounds))

  All_receptor = LR_data%>%
    dplyr::filter(!Receptor_genes=="")%>%
    dplyr::select(Compounds, Receptor_genes)
  All_receptor = split(All_receptor,list(All_receptor$Compounds))

  All_pre_receptor = LR_data%>%
    dplyr::filter(!precursor_transport=="")%>%
    dplyr::select(Compounds, precursor_transport)
  All_pre_receptor = split(All_pre_receptor,list(All_pre_receptor$Compounds))

  Synthetic_database = lapply(X = All_synthetic, FUN = function(x) {
    x <- unlist(strsplit(x$Synthetic_genes,split = ","))
    x <- x[nzchar(x)]
  })

  if (Synthetic_step2 == T){
  Synthetic_step2_database = lapply(X = All_synthetic_lv2, FUN = function(x) {
    x <- unlist(strsplit(x$Synthetic_genes_lv2,split = ","))
  })  } else {
    Synthetic_step2_database = NULL
  }

  Transport_database = lapply(X = All_trans, FUN = function(x) {
    x <- unlist(strsplit(x$Transporter_genes,split = ","))
  })

  Receptor_database = lapply(X = All_receptor, FUN = function(x) {
    x <- unlist(strsplit(x$Receptor_genes,split = ","))
  })

  if (precusor_transporter == T){
  Pre_trans_database = lapply(X = All_pre_receptor, FUN = function(x) {
    x <- unlist(strsplit(x$precursor_transport,split = ","))
  })  } else {
    Synthetic_step2_database = NULL
  }

  #Check genes with expression database
  allgenes = rownames(ave_expr)

  Synthetic_only = Synthetic_database[!names(Synthetic_database)%in%names(Transport_database)]
  Syn_both = Synthetic_database[names(Synthetic_database)%in%names(Transport_database)]
  Trans_both = Transport_database[names(Syn_both)]

  ### Expression of Only Syns
  syns_matrix = matrix(0, nrow = length(Synthetic_only), ncol = ncol(ave_expr))
  rownames(syns_matrix) = names(Synthetic_only)
  colnames(syns_matrix) = colnames(ave_expr)
  for (i in 1:length(Synthetic_only)) {
    synsg = Synthetic_only[[i]]
    Exp = matrix(NA, nrow = length(synsg), ncol = ncol(ave_expr))
    rownames(Exp) = synsg
    colnames(Exp) = colnames(ave_expr)
    for (j in 1:length(synsg)){
      gene = unlist(strsplit(synsg[[j]],split = " and "))
      gene_exp = ave_expr[rownames(ave_expr)%in%gene,,drop=F]
      Exp[j,] = apply(gene_exp, 2, function(x) exp(mean(log(x))))
    }
    Exp = na.omit(Exp)
    matrix = colMeans(Exp)
    if (is.null(Synthetic_step2_database)|names(Synthetic_only)[[i]] %notin% names(Synthetic_step2_database)){
      matrix_syng = matrix
    } else {
      synsg_step2 = c(unlist(Synthetic_step2_database[names(Synthetic_only)[[i]]]),unlist(Pre_trans_database[names(Synthetic_only)[[i]]]))
      Exp = matrix(NA, nrow = length(synsg_step2), ncol = ncol(ave_expr))
      rownames(Exp) = synsg_step2
      colnames(Exp) = colnames(ave_expr)
      for (j in 1:length(synsg_step2)){
        gene = unlist(strsplit(synsg_step2[[j]],split = " and "))
        gene_exp = ave_expr[rownames(ave_expr)%in%gene,,drop=F]
        Exp[j,] = apply(gene_exp, 2, function(x) exp(mean(log(x))))
      }
      Exp = na.omit(Exp)
      matrix_lv2 = colMeans(Exp)
      matrix_syng = sqrt(matrix * matrix_lv2)
    }
    syns_matrix[i,] = matrix_syng
  }

  ### Expression of Syns x Trans
  both_matrix = matrix(0, nrow = length(Syn_both), ncol = ncol(ave_expr))
  rownames(both_matrix) = names(Syn_both)
  colnames(both_matrix) = colnames(ave_expr)
  for (i in 1:length(Syn_both)) {
    #Synthetic
    synsg = Syn_both[[i]]
    Exp = matrix(NA, nrow = length(synsg), ncol = ncol(ave_expr))
    rownames(Exp) = synsg
    colnames(Exp) = colnames(ave_expr)
    for (j in 1:length(synsg)){
      gene = unlist(strsplit(synsg[[j]],split = " and "))
      gene_exp = ave_expr[rownames(ave_expr)%in%gene,,drop=F]
      Exp[j,] = apply(gene_exp, 2, function(x) exp(mean(log(x))))
    }
    Exp = na.omit(Exp)
    matrix = colMeans(Exp)
    if (is.null(Synthetic_step2_database)|names(Syn_both)[[i]] %notin% names(Synthetic_step2_database)){
      matrix_syng = matrix
    } else {
      synsg_step2 = c(unlist(Synthetic_step2_database[names(Syn_both)[[i]]]),unlist(Pre_trans_database[names(Syn_both)[[i]]]))
      Exp = matrix(NA, nrow = length(synsg_step2), ncol = ncol(ave_expr))
      rownames(Exp) = synsg_step2
      colnames(Exp) = colnames(ave_expr)
      for (j in 1:length(synsg_step2)){
        gene = unlist(strsplit(synsg_step2[[j]],split = " and "))
        gene_exp = ave_expr[rownames(ave_expr)%in%gene,,drop=F]
        Exp[j,] = apply(gene_exp, 2, function(x) exp(mean(log(x))))
      }
      Exp = na.omit(Exp)
      matrix_lv2 = colMeans(Exp)
      matrix_syng = sqrt(matrix * matrix_lv2)
    }

    transg = Trans_both[[i]]
    Exp = matrix(NA, nrow = length(transg), ncol = ncol(ave_expr))
    rownames(Exp) = transg
    colnames(Exp) = colnames(ave_expr)
    for (j in 1:length(transg)){
      gene = unlist(strsplit(transg[[j]],split = " and "))
      gene_exp = ave_expr[rownames(ave_expr)%in%gene,,drop=F]
      Exp[j,] = apply(gene_exp, 2, function(x) exp(mean(log(x))))
    }
    Exp = na.omit(Exp)
    matrix_trans = colMeans(Exp)

    both_matrix[i,] = sqrt(matrix_syng * matrix_trans)
  }

  all_matrix = rbind(syns_matrix,both_matrix)
  return(all_matrix)
}
