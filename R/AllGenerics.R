###########################
### annotation generics ###
###########################

#' @export
#' @rdname resetAnnotations
setGeneric("resetAnnotations", function(GIObject) standardGeneric("resetAnnotations"))

#' @export
#' @rdname annotateRegions
setGeneric("annotateRegions", function(GIObject, name, dat) standardGeneric("annotateRegions"))

#' @export
#' @rdname annotateInteractions
setGeneric("annotateInteractions",function(GIObject, annotations) standardGeneric("annotateInteractions"))

#######################
### export generics ###
#######################

#' @export
setGeneric("export.igraph",function(GIObject) standardGeneric("export.igraph"))

#' @export
setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts") standardGeneric("export.bed12"))

#' @export
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts") standardGeneric("export.chiasig"))

#######################
### getter generics ###
#######################

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

#######################
### helper generics ###
#######################

#' @rdname InteractionHelpers
#' @export
setGeneric("is.pp",function(GIObject) standardGeneric("is.pp"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.pd",function(GIObject) standardGeneric("is.pd"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.pt",function(GIObject) standardGeneric("is.pt"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.dd",function(GIObject) standardGeneric("is.dd"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.dt",function(GIObject) standardGeneric("is.dt"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.tt",function(GIObject) standardGeneric("is.tt"))

#' @rdname InteractionHelpers
#' @export
setGeneric("isInteractionType",function(GIObject, x, y) standardGeneric("isInteractionType"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.trans",function(GIObject) standardGeneric("is.trans"))

#' @rdname InteractionHelpers
#' @export
setGeneric("is.cis",function(GIObject) standardGeneric("is.cis"))

#######################
### setter generics ###
#######################

#' @rdname setters
#' @export
setGeneric("name<-",function(GIObject, value) standardGeneric("name<-"))

#' @rdname setters
#' @export
setGeneric("interactionCounts<-",function(GIObject, value) standardGeneric("interactionCounts<-"))

###########################
### processing generics ###
###########################

#' @export
setGeneric("countsBetweenAnchors",function(x, y, ...) standardGeneric ("countsBetweenAnchors"))

###########################
### subsetting generics ###
###########################

#' @export
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

##############################
### summarization generics ###
##############################

#' @export
setGeneric("summariseByFeatures", function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE)
    standardGeneric("summariseByFeatures"))

#' @export
setGeneric("summariseByFeaturePairs",function(GIObject, features.one, feature.name.one, features.two, feature.name.two)
    standardGeneric ("summariseByFeaturePairs"))

#############################
### miscellanous generics ###
#############################

#' @export
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE) standardGeneric("calculateDistances"))


