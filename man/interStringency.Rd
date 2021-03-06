% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interStringency.R
\name{interStringency}
\alias{interStringency}
\title{Extract interactions according to stringency and interaction type}
\usage{
interStringency(
  type = c("RBP", "ncRNA"),
  stringency = c("low", "medium", "high", "strict")
)
}
\arguments{
\item{type}{a character string indicating which interaction type is to be choosed,
. One of "RBP" (RBP-circRNA,RBP-lncRNA,miRNA-circRNA,miRNA-lncRNA,miRNA-RBP), "ncRNA (miRNA-circRNA,miRNA-lncRNA)": can be abbreviated}

\item{stringency}{a character string indicating which interaction stringency is to be choosed,
. One of "low" ( number of supported experiments > = 1  ), "medium ( > = 2)","high ( > = 3)","strict ( > = 5)"}
}
\value{
interaction of setting
}
\description{
interactions were extracted according to stringency and interaction type in the database of ENCORI
}
\examples{
\dontrun{
  data("dataM2C")
data("dataM2L")
data("dataM2R")
data("dataR2C")
data("dataR2L")
 interac <- interStringency(type = "ncRNA",stringency = "strict")
 interac <- interac[,c("node_gene_ID","target_gene_ID")]
}

}
