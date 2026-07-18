% Generated manually to match the current exported helper.
\name{Prepare_database}
\alias{Prepare_database}
\title{Create MetaLigand database lists}
\usage{
Prepare_database(LR_data = NULL, species = NULL)
}
\arguments{
\item{LR_data}{Optional ligand-receptor database table. If \code{NULL}, the
database bundled for \code{species} is loaded.}

\item{species}{Species database to use: \code{"mouse"}, \code{"human"}, or
\code{"zebrafish"}.}
}
\value{
A named list containing synthesis, transporter, receptor, and precursor
transporter gene databases.
}
\description{
Loads a MetaLigand ligand-receptor database and creates the database lists used
by \code{Meta_matrix()}.
}
