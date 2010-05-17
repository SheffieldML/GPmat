gpExtractParam <- function(model, only.values=TRUE) {

  if (model$learnScales) {
    fhandle <- get(paste(model$scaleTransform, "Transform", sep=""), mode="function")
    scaleParams = fhandle(model$scale, "xtoa")
    if (returnNames)
      for (i in 1:length(scaleParams))
	scaleParamNames[[i]] = paste("Output Scale ", as.character(i), sep="")
  } else {
    scaleParams = list()
    scaleParamNames = list()
  }

  ## Check if there is a mean function.
  meanFuncParams = list()
  if (("meanFunction" %in% names(model)) && !is.null(model$meanFunction))
    meanFuncParams = modelExtractParam(model$meanFunction, only.values)

  params <- kernExtractParam(model$kern, only.values)

  if (model.approx == "ftc") {
    params = c(kernParams, meanFuncParams, scaleParams)
    
    if (model$optimiseBeta) {
      fhandle <- get(paste(model$betaTransform, "Transform", sep=""), mode="function")
      betaParam = fhandle(model$beta, "xtoa")
      params = c(params, betaParam)
    }
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    paramPart = c(kernParams, meanFuncParams, scaleParams)
    if (model$optimiseBeta) {
      fhandle <- get(paste(model$betaTransform, "Transform", sep=""), mode="function")
      betaParam = fhandle(model$beta, "xtoa")
      paramPart = c(paramPart, betaParam)
    }
    if (model$fixInducing)
      params = paramPart
    else
      params = c(model$X_u, paramPart)
  }

  params <- Re(params)

  return (params)
}