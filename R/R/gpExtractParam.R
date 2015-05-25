gpExtractParam <- function(model, only.values=TRUE, ...) {

  ## Check if the output scales are being learnt.
  scaleParams = list()
  scaleParamNames = list()
  if (model$learnScales) {
    fhandle <- get(model$scaleTransform$func, mode="function")
    scaleParams = fhandle(model$scale, "xtoa")
    if (!only.values)
      for (i in 1:length(scaleParams))
	scaleParamNames[[i]] = paste("Output Scale ", as.character(i), sep="")
  }

  ## Check if there is a mean function.
  meanFuncParams = list()
  if (("meanFunction" %in% names(model)) && length(model$meanFunction)>0)
    meanFuncParams = modelExtractParam(model$meanFunction, only.values)

  kernParams <- kernExtractParam(model$kern, only.values)

  if (model$approx == "ftc") {
    params = unlist(c(kernParams, meanFuncParams, scaleParams))
    
    if (model$optimiseBeta) {
      fhandle <- get(model$betaTransform$func, mode="function")
      betaParam = fhandle(model$beta, "xtoa")
      params = c(params, betaParam)
    }
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    paramPart = unlist(c(kernParams, meanFuncParams, scaleParams))
    if (model$optimiseBeta) {
      fhandle <- get(model$betaTransform$func, mode="function")
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
