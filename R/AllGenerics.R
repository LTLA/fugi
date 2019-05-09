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
#' @rdname export.igraph
setGeneric("export.igraph",function(GIObject) standardGeneric("export.igraph"))

#' @export
#' @rdname export.bed12
setGeneric("export.bed12",function(GIObject, fn=NULL, score="counts") standardGeneric("export.bed12"))

#' @export
#' @rdname export.chiasig
setGeneric("export.chiasig", function(GIObject, fn=NULL, score="counts") standardGeneric("export.chiasig"))

#' @export
#' @rdname export.bedpe
setGeneric("export.bedpe", function(GIObject, fn=NULL, score="counts"){ standardGeneric("export.bedpe")} )

#######################
### getter generics ###
#######################

#' @export
#' @rdname getters
setGeneric("name", function(GIObject) standardGeneric("name"))

#' @export
#' @rdname getters
setGeneric("anchorOne",function(GIObject) standardGeneric("anchorOne"))

#' @export
#' @rdname getters
setGeneric("anchorTwo",function(GIObject) standardGeneric("anchorTwo"))

#' @export
#' @rdname getters
setGeneric("interactionCounts",function(GIObject) standardGeneric("interactionCounts"))

#' @export
#' @rdname getters
setGeneric("annotationFeatures",function(GIObject) standardGeneric("annotationFeatures"))

#######################
### helper generics ###
#######################

#' @export
#' @rdname InteractionHelpers
setGeneric("is.pp",function(GIObject) standardGeneric("is.pp"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.pd",function(GIObject) standardGeneric("is.pd"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.pt",function(GIObject) standardGeneric("is.pt"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.dd",function(GIObject) standardGeneric("is.dd"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.dt",function(GIObject) standardGeneric("is.dt"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.tt",function(GIObject) standardGeneric("is.tt"))

#' @export
#' @rdname InteractionHelpers
setGeneric("isInteractionType",function(GIObject, x, y) standardGeneric("isInteractionType"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.trans",function(GIObject) standardGeneric("is.trans"))

#' @export
#' @rdname InteractionHelpers
setGeneric("is.cis",function(GIObject) standardGeneric("is.cis"))

#######################
### setter generics ###
#######################

#' @export
#' @rdname setters
setGeneric("name<-",function(GIObject, value) standardGeneric("name<-"))

#' @export
#' @rdname setters
setGeneric("interactionCounts<-",function(GIObject, value) standardGeneric("interactionCounts<-"))

###########################
### processing generics ###
###########################

#' @export
#' @rdname countsBetweenAnchors-methods
setGeneric("countsBetweenAnchors",function(x, y, ...) standardGeneric ("countsBetweenAnchors"))

###########################
### subsetting generics ###
###########################

#' @export
#' @rdname subsetByFeatures
setGeneric("subsetByFeatures",function(GIObject, features, feature.class=NULL){standardGeneric ("subsetByFeatures")})

##############################
### summarization generics ###
##############################

#' @export
#' @rdname summariseByFeatures
setGeneric("summariseByFeatures", function(GIObject, features, feature.name, distance.method="midpoint", annotate.self=FALSE)
    standardGeneric("summariseByFeatures"))

#' @export
#' @rdname summariseByFeaturePairs
setGeneric("summariseByFeaturePairs",function(GIObject, features.one, feature.name.one, features.two, feature.name.two)
    standardGeneric ("summariseByFeaturePairs"))

#############################
### miscellanous generics ###
#############################

#' @export
#' @rdname calculateDistances
setGeneric("calculateDistances", function(GIObject, method="midpoint", floor=TRUE) standardGeneric("calculateDistances"))


