#' Functions to set interaction data 
#'
#' Use these functions to set data stored in the slots of a \linkS4class{GenomicInteractions} object.
#'
#' @name setters
#' @param GIObject A \linkS4class{GenomicInteractions} object
#' @param value A vector to replace a slot in the object
#' @return GenomicInteractions object
#' @rdname setters
#' @examples
#' anchor.one = GRanges(c("chr1", "chr1", "chr1", "chr1"), 
#'      IRanges(c(10, 20, 30, 20), width=5))
#' anchor.two = GRanges(c("chr1", "chr1", "chr1", "chr2"), 
#'      IRanges(c(100, 200, 300, 50), width=5))
#' interaction_counts = sample(1:10, 4)
#' test <- GenomicInteractions(anchor.one, anchor.two, 
#'      metadata=list(experiment_name="test", description="this is a test"),
#'      counts=interaction_counts)
#' metadata(test)
#'
#' name(test) <- "Mouse test"
#' name(test)
#'
#' description(test) <- "This is a test using the mouse genome"
#' description(test)
#'
#' interactionCounts(test) <- c(2,3,8,5)
#' interactionCounts(test)
NULL

#' @export
#' @rdname setters
#' @importFrom S4Vectors metadata<-
setReplaceMethod("name", "GenomicInteractions", function(GIObject, value){
    metadata(GIObject)$experiment_name = value
    GIObject
})

#' @rdname setters
#' @inheritParams Biobase::'description<-'
#' @importMethodsFrom Biobase 'description<-'
#' @importFrom S4Vectors metadata<-
#' @export
setReplaceMethod("description", "GenomicInteractions", function(object, value){
    metadata(object)$description = value
    object
})

#' @rdname setters
#' @export
#' @importFrom S4Vectors mcols<-
setReplaceMethod("interactionCounts", "GenomicInteractions", function(GIObject, value){
    if (!all(value == floor(value)))
        stop("value must contain integer values")
    value = as.integer(value)
    if (length(value) == 1)
        value = rep(value, length(GIObject))
    mcols(GIObject)$counts <- value
    GIObject
})


#' @export
#' @rdname setters
#' @inheritParams GenomeInfoDb::'seqinfo<-'
#' @importFrom GenomeInfoDb seqinfo<-
#' @importFrom GenomicInteractions regions regions<-
setMethod("seqinfo<-", "GenomicInteractions", function(x, new2old = NULL, 
    pruning.mode = c("error", "coarse", "fine", "tidy"), value) 
{
    regs <- regions(x, type="both")
    seqinfo(regs[[1]], new2old=new2old, pruning.mode=match.arg(pruning.mode)) <- value
    seqinfo(regs[[2]], new2old=new2old, pruning.mode=match.arg(pruning.mode)) <- value
    regions(x, type="both") <- regs
    x
})
