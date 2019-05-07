#' @title
#' Export interactions to an igraph object
#'
#' @description
#' Exports a \linkS4class{GenomicInteractions} object to graph object for use by \pkg{igraph} package. 
#' This uses unique anchors as nodes and generates edges between them. 
#'
#' For the resulting graph to be easily interpretable, anchors should be non-overlapping. 
#' This should already be the case for HiC data using either binned genomic regions or restriction fragments.
#' However, ChIA-PET data can contain overlapping anchors, which may need to be reduced to non-overlapping regions before graph export.
#' 
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#'
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' ig <- export.igraph(hic_example_data)
#' 
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @return A graph representation of the \linkS4class{GenomicInteractions} object.
#' @importFrom igraph graph_from_data_frame
#' @importFrom S4Vectors mcols
#' @importFrom GenomicInteractions anchors
#'
#' @export
#' @docType methods
#' @rdname export.igraph
#' @export
setMethod("export.igraph", "GenomicInteractions", function(GIObject){
    regs <- .get_single_regions(GIObject)
    nodes <- names(regs)
    if (is.null(nodes)) {
        nodes <- as.character(regs)
    }
    nodes <- data.frame(name = nodes, mcols(regs))

    a1 <- anchors(GIObject, type=1, id=TRUE)
    a2 <- anchors(GIObject, type=2, id=TRUE)
    edges <- data.frame(from = nodes$name[a1], to = nodes$name[a2])
    edges <- cbind(edges, mcols(GIObject))

    graph_from_data_frame(edges, directed=FALSE, vertices = nodes)
})
