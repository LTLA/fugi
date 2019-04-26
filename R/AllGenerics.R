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

