% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featurefinder.r
\name{getVarAv}
\alias{getVarAv}
\title{getVarAv}
\usage{
getVarAv(dd, varAv, varString)
}
\arguments{
\item{dd}{A dataframe}

\item{varAv}{A string corresponding to the numeric field to be averaged within each leaf node}

\item{varString}{A string}
}
\value{
An average of the numeric variable varString in the segment
}
\description{
This function generates a residual tree on a subset of the data
}
\examples{

require(featurefinder)
data(examples)
av=getVarAv(dat,"expected",pathterms)
}
\keyword{saveTree}
