#' @export
#' @rdname resetAnnotations
setGeneric("resetAnnotations", function(GIObject)standardGeneric ("resetAnnotations"))

#' @export
#' @rdname annotateRegions
setGeneric("annotateRegions", function(GIObject, name, dat) standardGeneric ("annotateRegions"))

#' @export
#' @rdname annotateInteractions
setGeneric("annotateInteractions",function(GIObject, annotations) standardGeneric("annotateInteractions"))

#' @export
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE) standardGeneric ("calculateDistances"))

#' @export
setGeneric("export.igraph",function(GIObject) standardGeneric ("export.igraph"))

#' @export
setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts") standardGeneric ("export.bed12"))

#' @export
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts") standardGeneric("export.chiasig"))

