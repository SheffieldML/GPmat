setClass("scoreList", 
	representation(params = "list", loglikelihoods = "numeric",
                       baseloglikelihoods = "numeric",
                       genes = "list", modelArgs = "list",
                       knownTargets = "character", TF = "character",
                       sharedModel = "list"),
         prototype(params=list(), loglikelihoods = numeric(),
                   baseloglikelihoods = numeric(),
                   genes=list(), modelArgs=list(), knownTargets="",
                   TF="", sharedModel=list()))


setMethod("initialize", "scoreList",
          function(.Object, params, loglikelihoods,
                   baseloglikelihoods=numeric(), genes, modelArgs,
                   knownTargets="", TF="", sharedModel=list()) {
            if (!is.list(sharedModel))
              sharedModel <- list(sharedModel)

            if (is.null(knownTargets))
              knownTargets <- ""

            if (is.null(TF))
              TF <- ""

            names(params) <- genes
            names(modelArgs) <- genes
            names(loglikelihoods) <- genes
            if (length(baseloglikelihoods) > 0)
              names(baseloglikelihoods) <- genes
  
            .Object@params <- params
            .Object@loglikelihoods <- loglikelihoods
            .Object@baseloglikelihoods <- baseloglikelihoods
            .Object@genes <- genes
            .Object@modelArgs <- modelArgs
            .Object@knownTargets <- knownTargets
            .Object@TF <- TF
            .Object@sharedModel <- sharedModel
            .Object
          })


setMethod("show", "scoreList",
          function(object) {
            if (length(object@loglikelihoods) != 1)
              genetext <- "genes"
            else
              genetext <- "gene"
            
            if (object@TF == "")
              cat("scoreList of ", length(object@loglikelihoods), " ", genetext, ".\n", sep="")
            else
              cat("scoreList of ", length(object@loglikelihoods), " ", genetext, " for TF ", object@TF, ".\n", sep="")
            if (all(object@knownTargets != "")) {
              cat("  Known targets: ", paste(object@knownTargets, collapse=", "), "\n")
            }
            idx <- .listSelectSomeIndex(object@genes, maxToShow=4)
            l <- object@genes[c(idx[[1]], idx[[3]]), drop=FALSE]
            itms <- c(l[idx[[1]]], idx[[2]],
                      if (!is.null(idx[[1]])) l[-idx[[1]]] else NULL)
            cat("  Genes: ", paste(itms, collapse=", "), sep="")
            if (length(object@genes)>length(itms))
              cat(" (",length(object@genes)," total)", sep="")
            cat("\n")
          })


setGeneric("loglikelihoods",    function(object) standardGeneric("loglikelihoods"))
setGeneric("loglikelihoods<-",  function(object, value) standardGeneric("loglikelihoods<-"))
setMethod("loglikelihoods", "scoreList", function(object) object@loglikelihoods)
setReplaceMethod("loglikelihoods", c("scoreList", "numeric"),
                 function(object, value) {
                   if (length(value) != length(object@loglikelihoods))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@loglikelihoods ), ")",sep=""))
                   object@loglikelihoods <- value
                   object
                 })

setGeneric("baseloglikelihoods",    function(object) standardGeneric("baseloglikelihoods"))
setGeneric("baseloglikelihoods<-",  function(object, value) standardGeneric("baseloglikelihoods<-"))
setMethod("baseloglikelihoods", "scoreList", function(object) object@baseloglikelihoods)
setReplaceMethod("baseloglikelihoods", c("scoreList", "numeric"),
                 function(object, value) {
                   if (length(value) != length(object@baseloglikelihoods))
                     stop(paste("the length of replacement (",
                                length(value),
                                ") should equal the existing length (",
                                length( object@baseloglikelihoods ), ")",sep=""))
                   object@baseloglikelihoods <- value
                   object
                 })

setGeneric("params",    function(object) standardGeneric("params"))
setGeneric("params<-",  function(object, value) standardGeneric("params<-"))
setMethod("params", "scoreList", function(object) object@params)
setReplaceMethod("params", c("scoreList", "list"),
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
setReplaceMethod("genes", c("scoreList", "list"),
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
setReplaceMethod("modelArgs", c("scoreList", "list"),
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
setReplaceMethod("knownTargets", c("scoreList", "character"),
                 function(object, value) {
                   object@knownTargets <- value
                   object
                 })

setGeneric("TF",    function(object) standardGeneric("TF"))
setGeneric("TF<-",  function(object, value) standardGeneric("TF<-"))
setMethod("TF", "scoreList", function(object) object@TF)
setReplaceMethod("TF", c("scoreList", "character"),
                 function(object, value) {
                   object@TF <- value
                   object
                 })

setGeneric("sharedModel",    function(object) standardGeneric("sharedModel"))
setGeneric("sharedModel<-",  function(object, value) standardGeneric("sharedModel<-"))
setMethod("sharedModel", "scoreList", function(object) object@sharedModel)
setReplaceMethod("sharedModel", c("scoreList", "list"),
                 function(object, value) {
                   object@sharedModel <- value
                   object
                 })

setMethod("length", "scoreList", function(x) length(x@loglikelihoods))

setMethod("[",
          signature(x="scoreList"),
          function(x, i, j, ..., drop) {
            par <- x@params[i]
            ll <- x@loglikelihoods[i]
            bll <- x@baseloglikelihoods[i]
            genes <- x@genes[i]
            args <- x@modelArgs[i]
            new("scoreList", params=par, loglikelihoods=ll,
                baseloglikelihoods=bll, genes=genes,
                modelArgs=args, knownTargets=x@knownTargets, TF=x@TF,
                sharedModel=x@sharedModel)
          })

setMethod("c", signature(x="scoreList"),
          function(x, ..., recursive=FALSE) {
            lists <- unlist(list(x, ...))
            params <- do.call(c, lapply(lists, function(y) y@params))
            loglikelihoods <- do.call(c, lapply(lists, function(y) y@loglikelihoods))
            baseloglikelihoods <- do.call(c, lapply(lists, function(y) y@baseloglikelihoods))
            genes <- do.call(c, lapply(lists, function(y) y@genes))
            modelArgs <- do.call(c, lapply(lists, function(y) y@modelArgs))
            new("scoreList", params=params, loglikelihoods=loglikelihoods,
                baseloglikelihoods=baseloglikelihoods, genes=genes,
                modelArgs=modelArgs, knownTargets=lists[[1]]@knownTargets,
                TF=lists[[1]]@TF, sharedModel=lists[[1]]@sharedModel)
          })

setGeneric("sort", function(x, decreasing=FALSE, ...) standardGeneric("sort"))
setMethod("sort", signature(x="scoreList"), 
          function(x, decreasing=FALSE, ...) {
            r <- sort(loglikelihoods(x), decreasing=decreasing, index.return=TRUE)
            x[r$ix]
          })

setGeneric("write.scores", function(x, ...) standardGeneric("write.scores"))
setMethod("write.scores", signature(x="scoreList"), 
          function(x, ...) {
            v <- cbind(loglikelihoods(x), baseloglikelihoods(x))
            colnames(v) <- c('log-likelihood', 'null_log-likelihood')
            write.table(v, ...)
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
setReplaceMethod("modelStruct", signature(object="GPModel",value="list"),
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
setReplaceMethod("var.exprs", "ExpressionTimeSeries",
                 function(object, value) assayDataElementReplace(object, "var.exprs", value))
