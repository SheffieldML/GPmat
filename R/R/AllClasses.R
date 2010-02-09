require(Biobase)

setClass("scoreList", 
	representation(params = "list", LLs = "numeric", genes = "list", modelArgs = "list", knownTargets = "character", TF = "character"))

	#prototype = list(params = list(), LLs = array(), genes = list(), useGPsim = array()))


scoreList <- function(params, LLs, genes, modelArgs, knownTargets, TF) {
  if (is.null(knownTargets))
    knownTargets <- ""

  if (is.null(TF))
    TF <- ""

  names(params) <- genes
  names(modelArgs) <- genes
  names(LLs) <- genes
  
  new("scoreList", params = params, LLs = LLs, genes = genes, modelArgs = modelArgs, knownTargets = knownTargets, TF = TF)
}


is.scoreList <- function(object) {
  return (class(object) == "scoreList")
}

setMethod("show", "scoreList",
          function(object) {
            if (object@TF == "")
              cat("Score list of", length(object@LLs), "genes.\n")
            else
              cat("Score list of", length(object@LLs), "genes for TF ", object@TF, ".\n")
            if (object@knownTargets != "") {
              cat("Known targets: ")
              print(object@knownTargets)
            }
          })

setClass("GPmodel", 
	representation(model = "list", type = "character"))


GPmodel <- function(model) {
  new("GPmodel", model = model, type = model$type)
}


is.GPmodel <- function(object) {
  return (class(object) == "GPmodel")
}

print.GPmodel <- function(m) {
  if (m@type == "cgpdisim")
    gpdisimDisplay(m@model$comp[[1]])
  else
    gpsimDisplay(m@model$comp[[1]])
  cat("  Log-likelihood:", modelLogLikelihood(m), "\n")
}

setMethod("show", "GPmodel",
          function(object) { print.GPmodel(object) })


setClass("GPdata", 
	representation(y = "list", yvar = "list", times = "array",
                       genes = "character", scale = "array", zScores = "array",
                       annotation = "character",
                       phenoData = "AnnotatedDataFrame",
                       featureData = "AnnotatedDataFrame"))


GPdata <- function(y, yvar, times, genes, scale, zScores, annotation, phenoData, featureData) {
  new("GPdata", y = y, yvar = yvar, times = times, genes = genes, scale = scale, zScores = zScores, annotation=annotation, phenoData=phenoData, featureData=featureData)
}

is.GPdata <- function(object) {
  return (class(object) == "GPdata")
}

print.GPdata <- function(x) {
  cat("An object of class GPdata:\n")
  show(x@featureData)
  show(x@phenoData)
  cat("Annotation: ", x@annotation)
}

setMethod("show", "GPdata", function(object) print.GPdata(object))
