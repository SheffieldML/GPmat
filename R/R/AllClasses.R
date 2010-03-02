setClass("scoreList", 
	representation(params = "list", LLs = "numeric",
                       genes = "list", modelArgs = "list",
                       knownTargets = "character", TF = "character",
                       sharedModel = "list"))

	#prototype = list(params = list(), LLs = array(), genes = list(), useGPsim = array()))


scoreList <- function(params, LLs, genes, modelArgs, knownTargets="", TF="",
                      sharedModel=list()) {
  if (is.null(knownTargets))
    knownTargets <- ""

  if (is.null(TF))
    TF <- ""

  if (!is.list(sharedModel))
    sharedModel <- list(sharedModel)

  names(params) <- genes
  names(modelArgs) <- genes
  names(LLs) <- genes
  
  new("scoreList", params = params, LLs = LLs, genes = genes,
      modelArgs = modelArgs, knownTargets = knownTargets, TF = TF,
      sharedModel = sharedModel)
}


setMethod("show", "scoreList",
          function(object) {
            if (length(object@LLs) != 1)
              genetext <- "genes"
            else
              genetext <- "gene"
            
            if (object@TF == "")
              cat("scoreList of ", length(object@LLs), " ", genetext, ".\n", sep="")
            else
              cat("scoreList of ", length(object@LLs), " ", genetext, " for TF ", object@TF, ".\n", sep="")
            if (all(object@knownTargets != "")) {
              cat("  Known targets: ", paste(object@knownTargets, collapse=", "), "\n")
            }
            idx <- listSelectSomeIndex(object@genes, maxToShow=4)
            l <- object@genes[c(idx[[1]], idx[[3]]), drop=FALSE]
            itms <- c(l[idx[[1]]], idx[[2]],
                      if (!is.null(idx[[1]])) l[-idx[[1]]] else NULL)
            cat("  Genes: ", paste(itms, collapse=", "), sep="")
            if (length(object@genes)>length(itms))
              cat(" (",length(object@genes)," total)", sep="")
          })


setGeneric("loglikelihoods",    function(object) standardGeneric("loglikelihoods"))
setGeneric("loglikelihoods<-",  function(object, value) standardGeneric("loglikelihoods<-"))
setMethod("loglikelihoods", "scoreList", function(object) object@LLs)
setReplaceMethod("loglikelihoods", c("scoreList", "ANY"),
                 function(object, value) {
                   if (length(value) != length(object@LLs))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@LLs ), ")",sep=""))
                   object@LLs <- value
                   object
                 })

setGeneric("params",    function(object) standardGeneric("params"))
setGeneric("params<-",  function(object, value) standardGeneric("params<-"))
setMethod("params", "scoreList", function(object) object@params)
setReplaceMethod("params", c("scoreList", "ANY"),
                 function(object, value) {
                   if (length(value) != length(object@params))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@params ), ")",sep=""))
                   object@params <- value
                   object
                 })

setGeneric("genes",    function(object) standardGeneric("genes"))
setGeneric("genes<-",  function(object, value) standardGeneric("genes<-"))
setMethod("genes", "scoreList", function(object) object@genes)
setReplaceMethod("genes", c("scoreList", "ANY"),
                 function(object, value) {
                   if (length(value) != length(object@genes))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@genes ), ")",sep=""))
                   object@genes <- value
                   object
                 })

setGeneric("modelArgs",    function(object) standardGeneric("modelArgs"))
setGeneric("modelArgs<-",  function(object, value) standardGeneric("modelArgs<-"))
setMethod("modelArgs", "scoreList", function(object) object@modelArgs)
setReplaceMethod("modelArgs", c("scoreList", "ANY"),
                 function(object, value) {
                   if (length(value) != length(object@modelArgs))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@modelArgs ), ")",sep=""))
                   object@modelArgs <- value
                   object
                 })

setGeneric("knownTargets",    function(object) standardGeneric("knownTargets"))
setGeneric("knownTargets<-",  function(object, value) standardGeneric("knownTargets<-"))
setMethod("knownTargets", "scoreList", function(object) object@knownTargets)
setReplaceMethod("knownTargets", c("scoreList", "ANY"),
                 function(object, value) {
                   object@knownTargets <- value
                   object
                 })

setGeneric("TF",    function(object) standardGeneric("TF"))
setGeneric("TF<-",  function(object, value) standardGeneric("TF<-"))
setMethod("TF", "scoreList", function(object) object@TF)
setReplaceMethod("TF", c("scoreList", "ANY"),
                 function(object, value) {
                   object@TF <- value
                   object
                 })

setGeneric("sharedModel",    function(object) standardGeneric("sharedModel"))
setGeneric("sharedModel<-",  function(object, value) standardGeneric("sharedModel<-"))
setMethod("sharedModel", "scoreList", function(object) object@sharedModel)
setReplaceMethod("sharedModel", c("scoreList", "ANY"),
                 function(object, value) {
                   object@sharedModel <- value
                   object
                 })

setMethod("length", "scoreList", function(x) length(x@LLs))

setMethod("[",
          signature(x="scoreList"),
          function(x, i, j, ..., drop) {
            par <- x@params[i]
            ll <- x@LLs[i]
            genes <- x@genes[i]
            args <- x@modelArgs[i]
            new("scoreList", params=par, LLs=ll, genes=genes,
                modelArgs=args, knownTargets=x@knownTargets, TF=x@TF,
                sharedModel=x@sharedModel)
          })

setMethod("c", signature(x="scoreList"),
          function(x, ..., recursive=FALSE) {
            lists <- unlist(list(x, ...))
            params <- do.call(c, lapply(lists, function(y) y@params))
            LLs <- do.call(c, lapply(lists, function(y) y@LLs))
            genes <- do.call(c, lapply(lists, function(y) y@genes))
            modelArgs <- do.call(c, lapply(lists, function(y) y@modelArgs))
            new("scoreList", params=params, LLs=LLs, genes=genes,
                modelArgs=modelArgs, knownTargets=lists[[1]]@knownTargets,
                TF=lists[[1]]@TF, sharedModel=lists[[1]]@sharedModel)
          })

setGeneric("sort", function(x, decreasing=FALSE, ...) standardGeneric("sort"))
setMethod("sort", signature(x="scoreList"), 
          function(x, decreasing=FALSE, ...) {
            r <- sort(x@LLs, decreasing, index.return=TRUE)
            x[r$ix]
          })

setClass("GPModel", 
         representation(model = "list", type = "character"),
         prototype(model=list(), type=""))

setMethod("initialize", "GPModel",
          function(.Object, model) {
            .Object@model <- model
            .Object@type <- model$type
            .Object
          })

is.GPModel <- function(object) {
  return (class(object) == "GPModel")
}

setMethod("show", "GPModel",
          function(object) {
            if (object@type == "cgpdisim")
              gpdisimDisplay(modelStruct(object)$comp[[1]])
            else
              gpsimDisplay(modelStruct(object)$comp[[1]])
            cat("  Log-likelihood:", modelLogLikelihood(object), "\n")
          })

setGeneric("modelStruct",    function(object) standardGeneric("modelStruct"))
setGeneric("modelStruct<-",  function(object, value) standardGeneric("modelStruct<-"))
setMethod("modelStruct", "GPModel", function(object) object@model)
setReplaceMethod("modelStruct", signature(object="GPModel",value="matrix"),
                 function(object, value) {
                   object@model <- value
                   object@type <- value$type
                   object
                 })

setGeneric("modelType",    function(object) standardGeneric("modelType"))
setMethod("modelType", "GPModel", function(object) object@model)


setClass("ExpressionTimeSeries", contains="ExpressionSet")

setMethod("initialize", "ExpressionTimeSeries", function(.Object, 
                                           assayData = assayDataNew(exprs = exprs,
                                             var.exprs = var.exprs),
                                           exprs = new("matrix"),
                                           var.exprs = new("matrix"), ...) { 
  if (!missing(assayData) && any(!missing(exprs), !missing(var.exprs))) {
    warning("using 'assayData'; ignoring 'exprs', 'var.exprs'") 
  }
  callNextMethod(.Object, assayData = assayData, ...) 
})

setValidity("ExpressionTimeSeries", function(object) {
  if ("experiments" %in% varLabels(object) && "modeltime" %in% varLabels(object))
    TRUE
  else
    "Missing phenoData annotation field(s)."
})

setGeneric("var.exprs",    function(object) standardGeneric("var.exprs"))
setGeneric("var.exprs<-",  function(object, value) standardGeneric("var.exprs<-"))
setMethod("var.exprs", "ExpressionTimeSeries", function(object) assayDataElement(object,"var.exprs"))
setReplaceMethod("var.exprs", signature(object="ExpressionTimeSeries",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "var.exprs", value))
