#' @title
#' Functions to access data in a \linkS4class{GenomicInteractions} object
#'
#' @description
#' These functions can be used to access data stored in the slots of a \linkS4class{GenomicInteractions} object.
#' They are mostly retained for convenience and backwards compatibility -
#' users are recommended to use the getters from the \pkg{GenomicInteractions} package itself.
#'
#' @name getters
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @rdname getters
#' 
#' @return 
#' For \code{anchorOne} and \code{anchorTwo}, a \linkS4class{GRanges} for the first and second anchor regions, respectively.
#'
#' For \code{interactionCounts}, a numeric vector with counts for each interaction in the object.
#'
#' For \code{description} and \code{name}, a string containing the description and name of the object, respectively.
#'
#' For \code{annotationFeatures}, a character vector of features with which the object was previously annotated, 
#' or \code{NA} if the object is unannotated.
#'
#' For \code{seqinfo}, a \linkS4class{Seqinfo} object for the regions in \code{GIObject}.
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @examples
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'      IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'      IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, 
#'      counts=interaction_counts)
#'
#' name(test)
#' description(test)
#' anchorOne(test)
#' anchorTwo(test)
#' interactionCounts(test)
#'
NULL

#' @rdname getters
#' @export
#' @aliases name
#' @importFrom S4Vectors metadata
setMethod("name", "GenomicInteractions", function(GIObject) metadata(GIObject)$experiment_name)

#' @rdname getters
#' @inheritParams Biobase::description
#' @importMethodsFrom Biobase description
#' @export
#' @importFrom S4Vectors metadata
setMethod("description", "GenomicInteractions", function(object) metadata(object)$description)

#' @rdname getters
#' @export
#' @importFrom GenomicInteractions anchors
setMethod("anchorOne", "GenomicInteractions", function(GIObject)  anchors(GIObject, type = 1)) 

#' @rdname getters
#' @export
#' @importFrom GenomicInteractions anchors
setMethod("anchorTwo", "GenomicInteractions", function(GIObject) anchors(GIObject, type = 2))

#' @rdname getters
#' @export
#' @importFrom S4Vectors mcols
setMethod("interactionCounts", "GenomicInteractions", function(GIObject){
    ## N.B. this may not apply to all GenomicInteractions objects
    if (!"counts" %in% names(mcols(GIObject))){
        warning("'counts' not in mcols of object; will return NULL")
    }
    mcols(GIObject)$counts
})

#' @importFrom GenomicInteractions regions
.get_single_regions <- function(GIObject) {
    regs <- regions(GIObject, type="both")
    if (!identical(regs[[1]], regs[[2]])) {
        stop("expecting the same set of regions for both anchors")
    }
    regs[[1]]
}

#' @rdname getters
#' @export
#' @importFrom S4Vectors mcols
setMethod("annotationFeatures", "GenomicInteractions", function(GIObject){
    regs <- .get_single_regions(GIObject)
    if ("node.class" %in% names(mcols(regs))) {
        annotation <- regs$node.class
    } else { 
        annotation <- NA_character_ 
    }
    annotation
})

.has_annotations <- function(x) {
    out <- annotationFeatures(x)
    length(out) > 1L || is.na(out)
}

#' @export
#' @rdname getters
#' @param x A \linkS4class{GenomicInteractions} object
#' @importFrom GenomeInfoDb seqinfo
setMethod("seqinfo", "GenomicInteractions", function(x) {
    regs <- .get_single_regions(x)
    seqinfo(regs)
})
