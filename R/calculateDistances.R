#' Calculate interaction distances
#'
#' This computes the linear distances between anchoring regions for each interaction,
#' simply by calling the \code{\link{pairdist}} function from \pkg{GenomicInteractions}.
#' 
#' @param GIObject A \linkS4class{GenomicInteractions} object
#' @param method String indicating how to calculate distances.
#' Should be \code{"midpoint"}, \code{"outer"} or \code{"inner"}. 
#' @param floor A logical specifying whether to round down distances to the nearest base pair. 
#' Defaults to \code{TRUE}.
#'
#' @return A numeric vector containing the distances between anchors/GRanges.
#' If \code{floor=TRUE}, this vector is integer.
#' Interactions on different chromosomes have distances set to \code{NA}.
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @seealso
#' \code{\link{pairdist}}, which this function calls.
#'
#' @examples
#' anchor.one <- GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'     IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two <- GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'     IRanges(c(100, 200, 300, 50), width=5))
#' test <- GenomicInteractions(anchor.one, anchor.two)
#' calculateDistances(test, method="midpoint")
#'         
#' @docType methods
#' @rdname calculateDistances
#' @export
#' @importFrom GenomicInteractions pairdist
setMethod("calculateDistances", "GenomicInteractions", function(GIObject, method="midpoint", floor=TRUE) { 
    if (method=="midpoint"){
        distances <- pairdist(GIObject, type = "mid")
    } else if(method=="outer"){
        distances <- pairdist(GIObject, type = "span")
    } else if(method=="inner"){
        distances <- pairdist(GIObject, type = "gap")
    }else{
        distances <- pairdist(GIObject, type = method)
    }
        
    if(floor) {
        floor(distances)
    } else{
        distances
    }
})
