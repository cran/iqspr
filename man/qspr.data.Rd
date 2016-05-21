\name{qspr.data}
\alias{qspr.data}
\docType{data}
\title{
sample dataset for QSPR modeling
}
\description{
This is sample dataset for QSPR modeling.
The first column is SMILES strings, and the rest 4 columns are properties. 
}
\usage{data("qspr.data")}
\examples{
data(qspr.data)
smis <- paste(qspr.data[,1])
ty <- qspr.data[,c(2,5)]
}
\keyword{datasets}
