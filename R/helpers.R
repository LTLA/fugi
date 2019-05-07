#' @title
#' Interaction type helper functions
#'
#' @description
#' Functions to classify interactions within \linkS4class{GenomicInteractions} objects.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' 
#' @details
#' \code{isInteractionType} identifies all interactions that occur between the annotated node classes \code{x} and \code{y}.
#'
#' \code{is.trans} and \code{is.cis} select trans-chromosomal and intra-chromosomal interactions, respectively.
#'
#' The other functions are convenience wrappers for identifying interactions between regions of common annotations,
#' namely promoter \code{p}, distal \code{d} or terminator \code{t} regions.
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @return A logical vector indicating whether each entry of \code{GIObject} is of the specified type.
#' @examples 
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' table(is.cis(hic_example_data))
#' sum(interactionCounts(hic_example_data))
#' 
#' @name InteractionHelpers
#' @rdname InteractionHelpers
NULL

#' @rdname InteractionHelpers
#' @export
setMethod("is.pp", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "promoter", "promoter"))

#' @rdname InteractionHelpers
#' @export
setMethod("is.pd", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "promoter", "distal"))

#' @rdname InteractionHelpers
#' @export
setMethod("is.pt", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "promoter", "terminator"))

#' @rdname InteractionHelpers
#' @export
setMethod("is.dd", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "distal", "distal"))

#' @rdname InteractionHelpers
#' @export
setMethod("is.dt", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "distal", "terminator"))

#' @rdname InteractionHelpers
#' @export
setMethod("is.tt", "GenomicInteractions", function(GIObject) isInteractionType(GIObject, "terminator", "terminal"))

#' @rdname InteractionHelpers
#' @param x,y Names of annotated node classes
#' @export
#' @importFrom GenomicInteractions anchors
setMethod("isInteractionType", "GenomicInteractions", function(GIObject, x, y) {
    anno <- annotationFeatures(GIObject)
    a1 <- anchors(GIObject, type=1, id=TRUE)
    a2 <- anchors(GIObject, type=2, id=TRUE)
    anno1 <- anno[a1]
    anno2 <- anno[a2]

    found <- anno1 %in% x & anno2 %in% y
    if (!identical(x, y)) {
        found <- found | (anno1 %in% y & anno2 %in% x)
    }
    found
})

#' @rdname InteractionHelpers
#' @export
setMethod("is.trans", "GenomicInteractions", function(GIObject) is.na(pairdist(GIObject)))

#' @rdname InteractionHelpers
#' @export
setMethod("is.cis", "GenomicInteractions", function(GIObject) !is.na(pairdist(GIObject)))

#' Return the total number of interactions in a GenomicInteractions GIObject
#'
#' @param x GenomicInteractions GIObject
#' @return The sum of the counts in GIObject
#' @docType methods
#' @export
setMethod("sum", "GenomicInteractions", function(x) sum(interactionCounts(x)))
