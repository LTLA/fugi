#' Reset annotations made to a GenomicInteractions object
#'
#' This function removes all annotations from a GenomicInteractions object by
#' deleting  all of the metadata columns associated with both anchors.
#'
#' @param GIObject An annotated GenomicInteractions object
#' @return invisible(1)
#' @examples
#' data(hic_example_data)
#' resetAnnotations(hic_example_data)
#'
#' @docType methods
#' @rdname resetAnnotations
#' @export
setGeneric("resetAnnotations", function(GIObject){standardGeneric ("resetAnnotations")})

#' @export
#' @rdname resetAnnotations
#' @importFrom S4Vectors mcols<- 
#' @importClassesFrom GenomicInteractions SimpleGenomicInteractions
setMethod("resetAnnotations", "SimpleGenomicInteractions", function(GIObject){ 
    objName <- deparse(substitute(GIObject))
    mcols(regions(GIObject)) <- NULL
    assign(objName, GIObject, envir = parent.frame())
    invisible(1)
})

#' Annotate regions
#'
#' Use this function to add metadata parallel to the `regions` slot of a 
#' GenomicInteractions or GenomicInteractions object.
#' 
#' @param GIObject A GenomicInteractions or GenomicInteractions object
#' @param name Character. Will be used as a column name.
#' @param dat Vector of the same length as the GenomicInteractions object,
#' containing data with which to annotate the object. 
#' @return invisible(1)
#' 
#' @rdname annotateRegions
#' @docType methods
#' @export
#' 
#' @examples 
#' data(hic_example_data)
#' chip <- runif(n = length(regions(hic_example_data)), max = 1000)
#' annotateRegions(hic_example_data, "chip", chip)
setGeneric("annotateRegions", function(GIObject, name, dat) standardGeneric ("annotateRegions"))

#' @export
#' @rdname annotateRegions
#' @importFrom S4Vectors mcols<- 
#' @importClassesFrom GenomicInteractions SimpleGenomicInteractions
setMethod("annotateRegions", c("GenomicInteractions", "character", "vector"), function(GIObject, name, dat) {
    objName <- deparse(substitute(GIObject))
    mcols(regions(GIObject)) <- dat
    assign(objName, GIObject, envir = parent.frame())
    invisible(1)
})

#' Annotate the interactions in a GenomicInteractions object
#'
#' This function will annotate both anchors with a list of named GRanges
#' objects. Each metadata column is labeled "name.id" and contains the id of
#' the genomic interval(s) it overlaps. Anonymous lists will be given names
#' "FEATURE#.id" where # is the position in the list.
#'
#' For each anchor a "node.class" metadata column will also be added, containing
#' the name of the list element which was \emph{first} annotated to each range.
#' Ranges with no overlaps will be classified as "distal". The identifiers for each 
#' individual feature/annotation are taken from either the name of the list item in the 
#' case of a GRangesList or from either the names of a the provided GRanges or an id column 
#' in its associated metadata.
#'
#' @param GIObject A GenomicInteractions object to be annotated
#' @param annotations A list containing GRanges (or GRangesList) objects with which to annotate
#'             the GenomicInteractions object.
#' @return invisible(1)
#' @rdname annotateInteractions
#' @docType methods
#' @export
#' 
#' @examples
#' 
#' library("GenomicRanges")
#' data(hic_example_data)
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl = split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
setGeneric("annotateInteractions",function(GIObject, annotations) standardGeneric("annotateInteractions"))

#' @export
#' @rdname annotateInteractions
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicInteractions GenomicInteractions
setMethod("annotateInteractions", c("GenomicInteractions", "list"), function(GIObject, annotations) {
    objName <- deparse(substitute(GIObject))
    mcols.reg <- mcols(regions(GIObject))
    mcols.reg$node.class <- NA
    if (is.null(names(annotations))){
        names(annotations) <- sprintf("FEATURE%i", seq_along(annotations))
    }

    feature_names_list <- lapply(annotations, .get_gr_names)
    if (any(vapply(feature_names_list, function(x) any(duplicated(x)), logical(1)))) {
        warning("Some features contain duplicate IDs which will result in duplicate annotations")
    }

    for(name in names(annotations)){
        message(paste("Annotating with", name, "..."))
        field_name <- paste(name, "id", sep=".")
        feature_names <- feature_names_list[[name]]
        mcols.reg[[field_name]] <- NA
        reg.ol <- findOverlaps(regions(GIObject), annotations[[name]])
        mcols.reg[[field_name]][ unique(queryHits(reg.ol)) ] <- split(feature_names[subjectHits(reg.ol)], 
                                                                      queryHits(reg.ol) )
        mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class) & !is.na(mcols.reg[[field_name]]), 
                                      name, mcols.reg$node.class)
    }

    mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class), 
                                  "distal", mcols.reg$node.class)
    mcols(regions(GIObject)) <- mcols.reg
    assign(objName, GIObject, envir = parent.frame())
    invisible(1)
})

#' @importClassesFrom GenomicRanges GenomicRanges GenomicRangesList
.get_gr_names <- function(x) {
    if (is(x, "GenomicRanges")) {
        if (!is.null(names(x))) {
            value <- names(x)
        } else if("id" %in% names(mcols(x))) {
            value <- x$id
        } else {
            stop("annotations requires an id column in elementMetadata or names to be non-null")
        }
    } else if(is(x, "GRangesList")) {
        value <- names(x)
    } else {
        stop("annotations must be GRanges or GRangesList objects")
    }
    as.character(value)
}
