#' @export
#' @rdname resetAnnotations
setGeneric("resetAnnotations", function(GIObject) standardGeneric("resetAnnotations"))

#' @export
#' @rdname annotateRegions
setGeneric("annotateRegions", function(GIObject, name, dat) standardGeneric("annotateRegions"))

#' @export
#' @rdname annotateInteractions
setGeneric("annotateInteractions",function(GIObject, annotations) standardGeneric("annotateInteractions"))

#' @export
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE) standardGeneric("calculateDistances"))

#' @export
setGeneric("export.igraph",function(GIObject) standardGeneric("export.igraph"))

#' @export
setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts") standardGeneric("export.bed12"))

#' @export
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts") standardGeneric("export.chiasig"))

#' @rdname getters
#' @export
setGeneric("name", function(GIObject) standardGeneric("name"))

#' @rdname getters
#' @export
setGeneric("anchorOne",function(GIObject) standardGeneric("anchorOne"))

#' @rdname getters
#' @export
setGeneric("anchorTwo",function(GIObject) standardGeneric("anchorTwo"))

#' @rdname getters
#' @export
setGeneric("interactionCounts",function(GIObject) standardGeneric("interactionCounts"))

#' @rdname getters
#' @export
setGeneric("annotationFeatures",function(GIObject) standardGeneric("annotationFeatures"))
