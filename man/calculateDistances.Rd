% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/calculateDistances.R
\docType{methods}
\name{calculateDistances}
\alias{calculateDistances}
\alias{calculateDistances,GenomicInteractions-method}
\title{Calculate interaction distances}
\usage{
calculateDistances(GIObject, method = "midpoint", floor = TRUE)

\S4method{calculateDistances}{GenomicInteractions}(GIObject,
  method = "midpoint", floor = TRUE)
}
\arguments{
\item{GIObject}{A \linkS4class{GenomicInteractions} object}

\item{method}{String indicating how to calculate distances.
Should be \code{"midpoint"}, \code{"outer"} or \code{"inner"}.}

\item{floor}{A logical specifying whether to round down distances to the nearest base pair. 
Defaults to \code{TRUE}.}
}
\value{
A numeric vector containing the distances between anchors/GRanges.
If \code{floor=TRUE}, this vector is integer.
Interactions on different chromosomes have distances set to \code{NA}.
}
\description{
This computes the linear distances between anchoring regions for each interaction,
simply by calling the \code{\link{pairdist}} function from \pkg{GenomicInteractions}.
}
\examples{
anchor.one <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
    IRanges(c(10, 20, 30, 20), width=5))
anchor.two <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
    IRanges(c(100, 200, 300, 50), width=5))
test <- GenomicInteractions(anchor.one, anchor.two)
calculateDistances(test, method="midpoint")
        
}
\seealso{
\code{\link{pairdist}}, which this function calls.
}
\author{
Malcolm Perry, Elizabeth Ing-Simmons
}
