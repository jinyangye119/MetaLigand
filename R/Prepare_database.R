#' Creat 4 databases Synthetic_database, Transport_database, Synthetic_step2_database, Pre_Receptor_database
#'
#' This function loads LR database csv file, and creat 4 basic database which used for "Mega_matrix" function
#'
#' @param infile "./raw"
#' @return A matrix of the infile
#' @export

Prepare_databse <- function(LR_data=NULL,
                            species = NULL){

  `%notin%` <- Negate(`%in%`)

  if (species%notin%c("mouse","human","zebrafish")){
  stop("species is required: mouse, human, zebrafish")
  }

  if (species == "mouse"){
    LR_data <- read.csv(system.file("extdata", "LR_data_mouse.csv", package = "MegaLigand"))
  } else if (species =="human"){
    LR_data <- read.csv(system.file("extdata", "LR_data_human.csv", package = "MegaLigand"))
  } else {
    LR_data <- read.csv(system.file("extdata", "LR_data_zebrafish.csv", package = "MegaLigand"))
  }

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
  Synthetic_step2_database = lapply(X = All_synthetic_lv2, FUN = function(x) {
    x <- unlist(strsplit(x$Synthetic_genes_lv2,split = ","))
  })
  Transport_database = lapply(X = All_trans, FUN = function(x) {
    x <- unlist(strsplit(x$Transporter_genes,split = ","))
  })
  Receptor_database = lapply(X = All_receptor, FUN = function(x) {
    x <- unlist(strsplit(x$Receptor_genes,split = ","))
  })
  precursor_transport_database = lapply(X = All_pre_receptor, FUN = function(x) {
    x <- unlist(strsplit(x$precursor_transport,split = ","))
  })

  All_database  = list(Synthetic_database,Synthetic_step2_database,Transport_database,Receptor_database,precursor_transport_database)
  names(All_database) = c("Synthetic_database","Synthetic_step2_database","Transport_database",
                          "Receptor_database","precursor_transport_database")

  return(All_database)
}

