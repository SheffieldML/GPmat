gpExpandParam <- function (model, params) {
  if (is.list(params)) {
    model$params = params
    params = params$xmin
  }

  if (model$approx == "ftc" || model$fixInducing)
    endVal = 0
  else {
    startVal = 1
    endVal = model$k * model$q
    model$X_u = matrix(params[startVal:endVal], model$k, model$q)
  }
  startVal = endVal + 1
  endVal = endVal + model$kern$nParams
  model$kern = kernExpandParam(model$kern, params[startVal:endVal])

  ## Check if there is a mean function.
  if ("meanFunction" %in% names(model) && length(model$meanFunction)>0) {
    startVal = endVal + 1
    endVal = endVal + model$meanFunction$numParams
    model$meanFunction = modelExpandParam(model$meanFunction, params[startVal:endVal])
  }

  ## Check if the output scales are being learnt.
  if (model$learnScales) {
    startVal = endVal + 1
    endVal = endVal + model$d
    fhandle <- get(model$scaleTransform$func, mode="function")
    model$scale = fhandle(params[startVal:endVal], "atox")
    model$m = gpComputeM(model)
  }

  ## Check if beta is being optimised.
  if (model$optimiseBeta) {
    startVal = endVal + 1
    endVal = endVal + prod(dim(as.matrix(model$beta)))
    fhandle <- get(model$betaTransform$func, mode="function")
    model$beta = fhandle(params[startVal:endVal], "atox")
  }

  ## Record the total number of parameters.
  model$nParams = endVal

  ## Update the kernel representations.
  if (model$approx == "ftc") {
    model = gpUpdateKernels(model, model$X, model$X_u)
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    model = gpUpdateKernels(model, model$X, model$X_u)
  } else
    stop("Unknown approximation type.")

  ## Update the vector 'alpha' for computing posterior mean.
  if ("alpha" %in% names(model))
    model = gpComputeAlpha(model)

  return (model)
}
