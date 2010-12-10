gpOut <- function(model, x) {
  if (nargs() < 2) {
    ## This implies evaluate for the training data.
    mu = model$mu ## bug: values of mu, varsigma not used
    varsigma = model$varSigma
  } else {
    if ('noise' %in% names(model)) {
      # [mu, varsigma] = gpPosteriorMeanVar(model, x)
      meanVar = gpPosteriorMeanVar(model, x, varsigma.return=TRUE)
      y = noiseOut(model$noise, meanVar$mu, meanVar$varsigma)
    } else
      y = gpPosteriorMeanVar(model, x)
  }

  return(y)
}
