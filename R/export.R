#' @title
#' Export interactions in BED12 format
#'
#' @description
#' Exports a \linkS4class{GenomicInteractions} object to file in BED12 format.
#' 
#' @details
#' BED12 files provide a method for visualising interactions.
#' However, it is not a good format for storing all of the data associated with an interaction dataset.
#' This is particularly true for trans-chromosomal interactions, which can only be stored in the BED12 \code{names} field.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param fn String containing a filename to write the object to.
#' @param score String or integer scalar specifying the metadata column to export as the \code{score} field.
#'
#' @return If \code{fn} is supplied, \code{NULL} is invisibly returned and a file is created at \code{fn}.
#'
#' Otherwise, a data.frame is returned containing all of the corresponding information from \code{GIObject}.
#' @export
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' export.bed12(hic_example_data, fn = tempfile(), score = "counts")
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @docType methods
#' @rdname export.bed12
#' @export
#' @importFrom utils write.table
#' @importFrom rtracklayer export
#' @importFrom grDevices col2rgb
setMethod("export.bed12", "GenomicInteractions", function(GIObject, fn=NULL, score="counts"){
    bed = asBED(GIObject, score)
    if (!is.null(fn)) {
        export(bed, fn, format="bed")
        return(invisible(NULL))
    } else {
        blocks = bed$blocks
        bed$blocks = NULL
        strand_info = as.character(strand(bed))
        strand_info[strand_info == "*"] = NA
        df = data.frame(
            seqnames=seqnames(bed),
            start=start(bed)-1,
            end=end(bed),
            names=bed$name,
            score=bed$score,
            strand=strand_info,
            thickStart=bed$thickStart,
            thickEnd=bed$thickEnd,
            itemRgb=apply(col2rgb(bed$itemRgb), 2, paste, collapse=","),
            blockCount=elementNROWS(blocks),
            blockStarts=unlist(lapply(start(blocks), paste, collapse = ","), use.names=FALSE),
            blockSizes=unlist(lapply(width(blocks), paste, collapse = ","), use.names=FALSE))
        return(bed)
    }
})

#' @title
#' Export interactions in BEDPE format.
#'
#' @description
#' Exports a \linkS4class{GenomicInteractions} object to file in the BEDPE (paired-end) format,.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param fn String containing a path to a file to write the interaction data to.
#' @param score String or integer scalar specifying the metadata column to export as the \code{score} field.
#'
#' @return If \code{fn} is supplied, \code{NULL} is invisibly returned and a file is created at \code{fn}.
#'
#' Otherwise, a data.frame is returned containing all of the corresponding information from \code{GIObject}.
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @export
#' @docType methods
#' @rdname export.bedpe
#' @export
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' export.bedpe(hic_example_data, fn = tempfile(), score = "counts")
setMethod("export.bedpe", "GenomicInteractions", function(GIObject, fn=NULL, score="counts"){
    score_vector = .getScore(GIObject, score)
    if (is.null(score_vector)) stop("Supplied score field not in element metadata.")
    output = cbind(as.character(seqnames(anchorOne(GIObject))),
                   start(anchorOne(GIObject))-1,
                   end(anchorOne(GIObject)),
                   as.character(seqnames(anchorTwo(GIObject))),
                   start(anchorTwo(GIObject))-1,
                   end(anchorTwo(GIObject)),
                   paste("interaction:", 1:length(GIObject), sep=""),
                   score_vector,
                   as.character(strand(anchorOne(GIObject))),
                   as.character(strand(anchorTwo(GIObject))))

    if(!is.null(fn)){
        write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
    }else{
        return(output)
    }

    return(invisible(1))
})

#' @title
#' Export interactions in a BEDPE-like format for use with ChiaSig
#'
#' @description
#' Exports a \linkS4class{GenomicInteractions} object to file in a BEDPE-like format for use with ChiaSig.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param fn String containing a path to a file to write the interaction data to.
#' @param score String or integer scalar specifying the metadata column to export as the \code{score} field.
#'
#' @return If \code{fn} is supplied, \code{NULL} is invisibly returned and a file is created at \code{fn}.
#'
#' Otherwise, a data.frame is returned containing all of the corresponding information from \code{GIObject}.
#'
#' @export
#' @docType methods
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' export.chiasig(hic_example_data, fn = tempfile(), score = "counts")
#'
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @rdname export.chiasig
#' @export
#' @references
#' Paulsen et al. (2014).
#' A statistical model of ChIA-PET data for accurate detection of chromatin 3D interactions.
#' \emph{Nucleic Acids Res.} 42(18), e143.
setMethod("export.chiasig", "GenomicInteractions", function(GIObject, fn=NULL, score="counts"){
    score_vec = .getScore(GIObject, score)
    if (is.null(score_vec)) stop("Supplied score field not in element metadata.")
    output = cbind(as.character(seqnames(anchorOne(GIObject))),
                    start(anchorOne(GIObject))-1,
                    end(anchorOne(GIObject)),
                    as.character(seqnames(anchorTwo(GIObject))),
                    start(anchorTwo(GIObject))-1,
                    score_vec)

    if (!is.null(fn)){
        write.table(output, fn, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE )
        return(invisible(NULL))
    } else{
        return(output)
    }
})

#' @importFrom S4Vectors mcols
.getScore <- function(x, score) {
    if (score=="counts") {
        ans <- interactionCounts(x)
    } else {
        ans <- mcols(x)[[score]]
    }
    if (is.null(ans)) {
        ans <- numeric(length(x))
    }
    ans
}

#' @importFrom S4Vectors mcols
.getNames <- function(x) {
    if ("name" %in% colnames(mcols(x))) {
        names <- mcols(x)[["name"]]
    } else {
        names <- sprintf("interaction_%i", seq_along(x))
    }
    names
}

#' @title Coerce to BED structure
#'
#' @description
#' Coerce the structure of a \linkS4class{GenomicInteractions} into a \linkS4class{GRanges} following BED-like conventions, 
#' i.e., with columns for blocks and thick regions.
#'
#' @param x A \linkS4class{GenomicInteractions} object.
#' @param keep.mcols Logical scalar, specifying whether to keep non-BED12 columns in the output object.
#' @param score String specifying which field to export as the \code{score} field in BED12.
#' Defaults to "auto" which will choose \code{"score"}, then \code{"counts"}, if present, or fill the column with zeros.
#'
#' @details
#' Note that setting \code{keep.mcols=TRUE} may cause problems with some parsers.
#'  
#' @return A \linkS4class{GRanges} object with the metadata columns \code{name}, \code{blockStarts} and \code{blockSizes}.
#'
#' @importFrom rtracklayer asBED
#' @importFrom IRanges PartitioningByWidth
#' @importFrom S4Vectors elementNROWS mcols
#' @importFrom BiocGenerics relist
#' @importFrom GenomicInteractions anchors swapAnchors
#' 
#' @author Malcolm Perry, Elizabeth Ing-Simmons
#' @rdname asBED
#' @export
#' @docType methods
#' @aliases asBED
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' asBED(hic_example_data)
setMethod("asBED", "GenomicInteractions", function(x, keep.mcols=FALSE, score="score") {
    if (!is.null(names(x))) 
        warning("Names will be dropped during BED12 export")

    x = swapAnchors(x, mode="order")
    a1 <- anchors(x, type=1, id=TRUE)
    a2 <- anchors(x, type=2, id=TRUE)

    # This is actually safe regardless of 'regions(x)'.
    # as swapAnchors standardizes both regions() anyway.
    # For the time being, at least; swapAnchors is not 
    # guaranteed to behave like this, so it might change.
    regs <- .get_single_regions(x) 

    is_trans = as.vector(seqnames(regs)[a1] != seqnames(regs)[a2])
    a1_cis = a1[!is_trans]
    a2_cis = a2[!is_trans]
    a1_trans = a1[is_trans]
    a2_trans = a2[is_trans]

    scores = .getScore(x, score)

    if (is.null(mcols(x)$color)) {
        mcols(x)$color = "#000000"
    }

    names = paste0(seqnames(regs)[a1], ":", start(regs)[a1] - 1 , "..", end(regs)[a1], "-",
        seqnames(regs)[a2], ":", start(regs)[a2] - 1, "..", end(regs)[a2], ",", scores)

    cis_blocks = relist(IRanges(
            start=c(rbind(rep(1L, length(a1_cis)), start(regs)[a2_cis] - start(regs)[a1_cis] + 1L)),
            width=c(rbind(width(regs)[a1_cis], width(regs)[a2_cis]))),
        PartitioningByWidth(rep(2, length(a1_cis))))

    output_cis = GRanges(
        seqnames=as.character(seqnames(regs)[a1_cis]),
        IRanges(start=start(regs)[a1_cis],
                end=end(regs)[a2_cis]),
        name=names[!is_trans],
        score=scores[!is_trans],
        strand=ifelse(
            strand(regs)[a1_cis] == strand(regs)[a2_cis] &
                   as.vector(strand(regs)[a1_cis]) %in% c("+", "-"),
            as.vector(strand(regs)[a1_cis]),
            "*"),
        thickStart=start(regs)[a1_cis],
        thickEnd=end(regs)[a2_cis],
        itemRgb=mcols(x)$color[!is_trans],
        blocks=cis_blocks
    )

    trans_blocks = relist(IRanges(
            start=rep(1, 2*length(a1_trans)),
            width=c(width(regs)[a1_trans], width(regs)[a2_trans])),
        PartitioningByWidth(rep(1, 2*length(a1_trans))))

    output_trans = GRanges(
        seqnames=c(as.character(seqnames(regs)[a1_trans]),
                as.character(seqnames(regs)[a2_trans])),
        IRanges(start=c(start(regs)[a1_trans],
                        start(regs)[a2_trans]),
                end=c(end(regs)[a1_trans],
                        end(regs)[a2_trans])),
        name=rep(names[is_trans], 2),
        score=rep(scores[is_trans], 2),
        strand=c(as.character(strand(regs)[a1_trans]),
                    as.character(strand(regs)[a2_trans])),
        thickStart=c(start(regs)[a1_trans],
                     start(regs)[a2_trans]),
        thickEnd=c(end(regs)[a1_trans],
                   end(regs)[a2_trans]),
        itemRgb=rep(mcols(x)$color[is_trans], 2),
        blocks=trans_blocks
    )

    extra_cols = setdiff(colnames(mcols(x)), c("score", "name"))

    if(length(extra_cols) && keep.mcols==TRUE) {
        mcols(output_cis) = cbind(mcols(output_cis), mcols(x)[!is_trans, extra_cols,drop=FALSE])
        mcols(output_trans) =
            cbind(mcols(output_trans),
                  rep(mcols(x)[is_trans, extra_cols,drop=FALSE], 2))
    }

    sort(c(output_cis, output_trans))
})
