#' Reset annotations
#'
#' This function removes all annotations from a \linkS4class{GenomicInteractions} object,
#' by deleting all of the metadata columns associated with both anchors.
#'
#' @param GIObject An annotated \linkS4class{GenomicInteractions} object.
#'
#' @return \code{NULL} is invisibly returned, 
#' while the value passed as \code{GIObject} is modified in the global namespace.
#' 
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' mcols(regions(hic_example_data, type=1))
#'
#' resetAnnotations(hic_example_data)
#' mcols(regions(hic_example_data, type=1))
#'
#' @docType methods
#' @export
#' @rdname resetAnnotations
#' @importFrom S4Vectors mcols<- 
#' @importClassesFrom GenomicInteractions GenomicInteractions 
#' @importFrom GenomicInteractions regions regions<-
setMethod("resetAnnotations", "GenomicInteractions", function(GIObject){ 
    objName <- deparse(substitute(GIObject))
    fsets <- regions(GIObject, type="both") 
    for (i in seq_along(fsets)) {
        mcols(fsets[[i]]) <- NULL
    }
    regions(GIObject, type="both") <- fsets
    assign(objName, GIObject, envir = parent.frame())
    invisible(NULL)
})

#' Annotate regions
#'
#' This function adds metadata parallel to the \code{regions} slot of a \linkS4class{GenomicInteractions} object.
#' 
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param name Character vector to be used as the column name.
#' @param dat Vector of the same length as \code{regions(GIObject, 1)},
#' containing data with which to annotate the object. 
#'
#' @details
#' This function only works if both \code{regions(GIObject, 1)} and \code{regions(GIObject, 2)} are identical.
#' In cases where they are not identical, it doesn't make sense to assign a single \code{dat} vector to the metadata.
#' 
#' @return \code{NULL} is invisibly returned, 
#' while the value passed as \code{GIObject} is modified in the global namespace.
#' 
#' @docType methods
#' @examples 
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' mcols(regions(hic_example_data, type=1))
#'
#' chip <- runif(length(regions(hic_example_data, type=1)), max=1000)
#' annotateRegions(hic_example_data, "chip", chip)
#' mcols(regions(hic_example_data, type=1))
#'
#' @export
#' @rdname annotateRegions
#' @importFrom S4Vectors mcols<-  List
#' @importClassesFrom GenomicInteractions GenomicInteractions 
#' @importFrom GenomicInteractions regions regions<-
setMethod("annotateRegions", c("GenomicInteractions", "character", "vector"), function(GIObject, name, dat) {
    objName <- deparse(substitute(GIObject))

    regs <- .get_single_regions(GIObject)
    mcols(regs)[[name]] <- dat
    regions(GIObject, type="both") <- List(regs, regs)

    assign(objName, GIObject, envir = parent.frame())
    invisible(1)
})

#' Annotate interactions 
#'
#' This function will annotate the interactions in a \linkS4class{GenomicInteractions} object with node class information.
#' 
#' @details
#' For each interaction, this function identifies the entries of \code{annotations} that overlap an anchor region.
#' For each entry in \code{annotations}, a metadata column is added to the \linkS4class{GenomicInteractions} object,
#' containing a list that specifies the elements inside the entry that overlap each interaction.
#' The metadata column is named as \code{NAME.id} where \code{NAME} is the name of the entry.
#' Anonymous lists will be given names \code{FEATURE#.id} where \code{#} is the position in the list.
#'
#' For each anchor, a \code{node.class} metadata column will also be added, 
#' containing the name of the list element which was \emph{first} annotated to each range.
#' Ranges with no overlaps will be classified as \code{"distal"}. 
#' The identifiers for each individual feature/annotation are taken from either the name of the list item in the case of a GRangesList or from either the names of a the provided GRanges or an id column in its associated metadata.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object to be annotated.
#' @param annotations A list containing \linkS4class{GRanges} (or \linkS4class{GRangesList}) objects to be used to annotate \code{GIObject}.
#'
#' @return \code{NULL} is invisibly returned, 
#' while the value passed as \code{GIObject} is modified in the global namespace.
#'
#' @docType methods
#' 
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' mcols(regions(hic_example_data, type=1))
#'
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl <- split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
#' mcols(regions(hic_example_data, type=1))
#'
#' @export
#' @rdname annotateInteractions
#' @importFrom S4Vectors queryHits subjectHits List
#' @importFrom GenomicInteractions GenomicInteractions 
#' @importFrom GenomicInteractions regions regions<-
setMethod("annotateInteractions", c("GenomicInteractions", "list"), function(GIObject, annotations) {
    objName <- deparse(substitute(GIObject))

    regs <- .get_single_regions(GIObject)

    mcols.reg <- mcols(regs)
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
        reg.ol <- findOverlaps(regs, annotations[[name]])
        mcols.reg[[field_name]][ unique(queryHits(reg.ol)) ] <- split(feature_names[subjectHits(reg.ol)], 
                                                                      queryHits(reg.ol) )
        mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class) & !is.na(mcols.reg[[field_name]]), 
                                      name, mcols.reg$node.class)
    }

    mcols.reg$node.class <- ifelse(is.na(mcols.reg$node.class), 
                                  "distal", mcols.reg$node.class)
    mcols(regs) <- mcols.reg

    regions(GIObject, type="both") <- List(regs, regs)
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
