#' Plot summary statistics 
#'
#' Makes summary plots of the counts, interaction distances, interaction annotations, 
#' and percentage of cis and trans interactions for a \linkS4class{GenomicInteractions} object.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param other Numeric scalar passed to \code{\link{plotInteractionAnnotations}}. 
#' @param cut Numeric scalar passed to \code{\link{plotCounts}}.
#'
#' @return A plot is made on the current graphics device.
#' The output of \code{\link{marrangeGrob}} is returned. 
#' @import ggplot2
#' @importFrom gridExtra marrangeGrob
#' @export
#'
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' plotSummaryStats(hic_example_data)
plotSummaryStats <- function(GIObject, other=5, cut=10){
  p1 <- plotCisTrans(GIObject)
  p2 <- plotDists(GIObject)
  
  if(!is.null(interactionCounts(GIObject))){
    p4 <- plotCounts(GIObject, cut=cut)
  } else {
    p4 <- NULL
  }
  
  if (.has_annotations(GIObject)) {
    p3 <- plotInteractionAnnotations(GIObject, other=other)
  }else{
    p3 <- NULL
  }
  
  p_list <- list(p4, p2, p1, p3)
  p_list <- p_list[lengths(p_list)>0]
  marrangeGrob(p_list, ncol=2, nrow =2, top = name(GIObject))
}

#' Plot cis/trans percentages
#'
#' Plots the percentages of cis and trans interactions for a \linkS4class{GenomicInteractions} object as a donut plot.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#'
#' @return A ggplot2 plot.
#' @export
#'
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' plotCisTrans(hic_example_data)
plotCisTrans <- function(GIObject){

  cis_p <- sum(is.cis(GIObject))
  trans_p <- sum(is.trans(GIObject))

  if (cis_p + trans_p != length(GIObject)){
    stop("cis + trans does not equal total interactions")
  }

  dat <- data.frame(count=c(cis_p, trans_p), category=c("cis", "trans"))
  dat$fraction = dat$count / sum(dat$count)
  dat = dat[order(dat$fraction), ]
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  dat$label <- paste(dat$category, "\n", signif(100*dat$fraction,3), "%")
  dat$labelpos <- (dat$ymax+dat$ymin)/2
  
  p <- ggplot(dat, aes_string(fill="category", ymax="ymax", ymin="ymin", 
                              xmax=4, xmin=2.3, label="label")) +
    geom_rect() +
    geom_text(aes_string(x=(4+2.3)/2, y="labelpos")) +
    coord_polar(theta="y") +
    xlim(c(0, 4))+
    theme_bw(base_size=16) +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(axis.title=element_blank()) +
    theme(panel.border=element_blank()) +
    theme(legend.position="none") +
    ggtitle("Cis / Trans Percentages")

  return(p)

}

#' Cis interaction distances
#'
#' Plots a histogram of interaction distances for a \linkS4class{GenomicInteractions} object.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param breaks A numeric vector of breaks for the histogram.
#' @param method String specifying the method used for distance between anchors, see \code{\link{calculateDistances}}.
#'
#' @return A ggplot2 plot.
#' @export
#'
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#' plotDists(hic_example_data)
plotDists <- function(GIObject, breaks=c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 2000000), method="midpoint"){

  options("scipen"=100, "digits"=4)
  dists <- calculateDistances(GIObject, method=method)
  dists <- dists[!is.na(dists)]

  breaks <- unique(sort(c(breaks, max(dists))))

  labs <- vector()
  for(i in 1:length(breaks)-2){
    labs[i] <- paste(breaks[i], breaks[i+1], sep="-")
  }
  labs[length(breaks)-1] <- paste(">", breaks[length(breaks)-1])

  dists_df <- data.frame(dists, Distance=cut(dists, breaks, labels=labs))
  p <- ggplot(dists_df, aes_string(x="Distance")) +
    geom_bar() +
    theme_bw(base_size=16)+
    theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
    xlab("Interaction distance (bp)") + ylab("Count")

  return(p)
}

#' Plot interaction types
#'
#' Plot a donut plot of interaction types for an annotated \linkS4class{GenomicInteractions} object.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object
#' @param node.classes Character vector specifying the node classes to include in the plot.
#' Defaults to all available node classes.
#' @param viewpoints Character vector of acceptable node classes,
#' where interactions are only considered if least one anchor is of an acceptable node class. 
#' Defaults to all classes in \code{node.classes}.
#' @param other Numeric scalar of the minimum proportion for an interaction type to be shown.
#' Interaction types making up fewer than \code{other} percent of the total interactions will be consolidated into a single \dQuote{other} category.
#' @param keep.order Logical scalar specifying whether the original order of \code{node.classes} should be retained for plotting.
#' Default is to use the alphabetical order.
#' @param legend Logical scalar specifying whether the legend should be plotted.
#' If \code{FALSE}, donut plot is annotated with category names.
#'
#' @return A ggplot2 plot.
#'
#' @export
#'
#' @examples
#' library("GenomicRanges")
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#'
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl = split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
#' plotInteractionAnnotations(hic_example_data)
plotInteractionAnnotations <- function(GIObject, node.classes=NULL, viewpoints=NULL, other=0, keep.order=FALSE, legend=FALSE){
  #add check for node classes existing!!
  dat <- categoriseInteractions(GIObject, node.classes, viewpoints)

  dat$fraction = dat$count / sum(dat$count)
  if(keep.order==FALSE) {
    dat = dat[order(dat$fraction), ]
  } else {
    dat = dat[nrow(dat):1, ] # order of node.classes specified in arguments
  }

  if (other > 0) {
    count.other = sum(dat$count[dat$fraction < other/100])
    fraction.other = sum(dat$fraction[dat$fraction < other/100])
    dat = dat[(dat$fraction > other/100),]
    rownames(dat) = 1:nrow(dat)
    other.row = data.frame("other", count.other, fraction.other)
    colnames(other.row) = colnames(dat)
    dat = rbind(other.row, dat) # needs to be at top to be plotted 'last' (anticlockwise)
  }

  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  dat$label <- paste(dat$category, "\n", signif(100*dat$fraction,3), "%")
  dat[dat$fraction < 0.03, "label"] <- NA
  dat$labelpos <- (dat$ymax+dat$ymin)/2
  dat$fraclabel <- paste(signif(100*dat$fraction,3), "%")
  
  p <- ggplot(dat, aes_string(fill="category", ymax="ymax", ymin="ymin", 
                              xmax=4, xmin=2.3)) +
    geom_rect() +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    ylim(c(0,1)) +
    theme_bw(base_size=16) +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(axis.title=element_blank()) +
    theme(panel.border=element_blank()) +
    ggtitle("Interaction Classes")

  if (!legend){
    p <- p + geom_text(aes_string(x=(4+2.3)/2, y="labelpos", label="label")) + 
      theme(legend.position="none")
  }else{
    p <- p + geom_text(aes_string(x=(4+2.3)/2, y="labelpos", label="fraclabel"))
  }

  return(p)
}

#' Categorize interaction types
#'
#' Get the numbers of interaction types existing in a \linkS4class{GenomicInteractions} object.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param node.classes Character vector specifying the node classes to include in the plot.
#' Defaults to all node classes.
#' @param viewpoints Character vector of acceptable node classes,
#' where interactions are only considered if least one anchor is of an acceptable node class. 
#' Defaults to all classes in \code{node.classes}.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
#' library("GenomicRanges")
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#'
#' data(mm9_refseq_promoters)
#' mm9_refseq_grl = split(mm9_refseq_promoters, mm9_refseq_promoters$id)
#' annotateInteractions(hic_example_data, list(promoter=mm9_refseq_grl))
#' categoriseInteractions(hic_example_data)
categoriseInteractions <- function(GIObject, node.classes=NULL, viewpoints=NULL){
  if(!.has_annotations(GIObject)) {
    stop("object not annotated!")
  }
  
  if(is.null(node.classes)){
    node.classes <- unique(annotationFeatures(GIObject))
  }

  if(is.null(viewpoints)) {
    viewpoints = node.classes
  } else if (!all(viewpoints %in% node.classes)) {
    stop("These viewpoints are not present in the dataset given")
  }

  line_count <- sum((length(node.classes):1)[1:length(viewpoints)])
  dat <- data.frame(category = rep(0, line_count), count=rep(0, line_count))

  line <- 0
  for (x in 1:length(viewpoints)){
    for (y in x:length(node.classes)){
      line <- line + 1
      dat[line,"category"] <- paste(viewpoints[x],node.classes[y], sep="-")
      dat[line, "count"] <- sum(isInteractionType(GIObject, viewpoints[x], node.classes[y]))
    }
  }
  return(dat)
}

#' Plot interaction count frequency
#'
#' Plot a bar chart of the number of interactions supported by different numbers of reads in a \linkS4class{GenomicInteractions} object.
#'
#' @param GIObject A \linkS4class{GenomicInteractions} object.
#' @param normalise Logical scalar specifying whether to plot the proportion of total reads instead of the count.
#' @param cut Numeric scalar specifying the upper limit on the count,
#' above which all interactions are aggregated into a single category.
#' Setting this to \code{NULL} will not perform any aggregation.
#'
#' @return A ggplot2 plot.
#' @export
#'
#' @examples
#' data(hic_example_data)
#' hic_example_data <- updateObject(hic_example_data)
#'
#' plotCounts(hic_example_data)
#' plotCounts(hic_example_data, normalise=TRUE)
plotCounts <- function(GIObject, normalise=FALSE, cut = 10){
  counts <- interactionCounts(GIObject)
  if(is.null(counts)){ stop("no count data to plot")}
  dat <- as.data.frame(table(counts), stringsAsFactors = FALSE)
  ylabel <- "Count"

  if (normalise){
    dat$Freq <- dat$Freq / sum(dat$Freq)
    ylabel <- "Proportion of all interactions"
  }

  if (!is.null(cut)){
    above_cut_sum <- sum(dat$Freq[as.numeric(dat$counts) > cut])
    dat <- dat[(as.numeric(dat$counts) <= cut),]
    dat <- rbind(dat, c(paste(">", cut), above_cut_sum))
    dat$Freq <- as.numeric(dat$Freq)
    dat$counts <- factor(dat$counts, levels=c(as.numeric(dat$counts[1:nrow(dat)-1]), paste(">", cut)))
  }else{
    dat$counts <- factor(dat$counts, levels=c(as.numeric(dat$counts)))
  }

  p <- ggplot(dat, aes_string(x="counts", y="Freq")) +
    geom_bar(stat="identity", position="dodge") +
    xlab("Number of reads") +
    ylab(ylabel) +
    theme_bw(base_size=16)
  return(p)
}
