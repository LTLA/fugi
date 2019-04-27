#' Functions to access data held in a GenomicInteractions object.
#'
#' Use these functions to access data stored in each of the slots of a
#' GenomicInteractions object.
#'
#' @name getters
#' @param GIObject A Gnteractions object
#' @rdname getters
#' 
#' @return For 'anchorOne' and 'anchorTwo', a GRanges. For 'interactionCounts', 
#' a numeric vector with counts for each interaction in the object. For
#'   'description' and 'name',  a character vector with
#'   length 1. For 'annotationFeatures', a character vector of features with
#'   which the object was previously annotated, or 'NA' if the object is unannotated.
#'
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

#' @rdname getters
#' @export
#' @importFrom IndexedRelations featureSets
#' @importFrom S4Vectors mcols
setMethod("annotationFeatures", "GenomicInteractions", function(GIObject){
    regs <- featureSets(GIObject)
    if (length(regs)==2L) {
        stop("expecting one set of regions only")
    }

    if ("node.class" %in% names(mcols(regs))) {
        annotation <- regs$node.class
    } else { 
        annotation <- NA_character_ 
    }

    annotation
})
