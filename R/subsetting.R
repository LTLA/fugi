#' Subset a \linkS4class{GenomicInteractions} object by features
#'
#' Subsets interactions for which at least one of the anchors overlaps with a given GRanges object.
#' Alternatively, subsets interactions based on annotated feature IDs for a particular feature.
#'
#' @docType methods
#' @param GIObject A \linkS4class{GenomicInteractions} object
#' @param features A GRanges or GRangesList object, or a character vector containing
#'                      IDs of annotated features, e.g. promoter IDs.
#' @param feature.class If `features' is a character vector, the corresponding feature name, e.g. "promoter".
#' @return a subsetted \linkS4class{GenomicInteractions} object
#' 
#' @examples 
#' data("hic_example_data")
#' hic_example_data <- updateObject(hic_example_data)
#'
#' data("mm9_refseq_promoters")
#' annotateInteractions(hic_example_data, list(promoter = mm9_refseq_promoters))
#' ids <- names(mm9_refseq_promoters[1:10])
#' subsetByFeatures(hic_example_data, ids, "promoter")
#'
#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @export
#' @importFrom IRanges overlapsAny
setMethod("subsetByFeatures", c("GenomicInteractions", "GRanges", "missing"), function(GIObject, features, feature.class=NULL){
    i <- overlapsAny(GIObject, features)
    GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @export
#' @importFrom IRanges overlapsAny
setMethod("subsetByFeatures", c("GenomicInteractions", "GRangesList", "missing"), function(GIObject, features, feature.class=NULL){
  i <- overlapsAny(GIObject, features)
  GIObject[i]
})

#' @rdname GenomicInteractions-subsetByFeatures-methods
#' @export
#' @importFrom GenomicInteractions anchors regions
#' @importFrom S4Vectors mcols
setMethod("subsetByFeatures", c("GenomicInteractions", "character", "character"), 
          function(GIObject, features, feature.class){
    if(!.has_annotations(GIObject) || !feature.class %in% annotationFeatures(GIObject)) {
      stop(paste(feature.class," has not been annotated on this object"))
    }

    #get regions which are annotated with given feature IDs
    FUN <- function(GIObject, type) {
        regs <- regions(GIObject, type=type)
        ids <- mcols(regs)[[paste(feature.class, "id", sep=".")]]
        which(sapply(ids, function(x) any(features %in% x)))
    }
    idx1 <- FUN(GIObject, 1)
    idx2 <- FUN(GIObject, 2)

    #get object index for region idx
    a1 <- anchors(GIObject, type=1, id=TRUE)
    a2 <- anchors(GIObject, type=2, id=TRUE)
    gi_idx <- (a1 %in% idx1) | (a2 %in% idx2)
    
    GIObject[gi_idx]
})
