gpComputeM <- function(model) {
  ## Remove mean function value from m (if mean function present).
  if ("meanFunction" %in% names(model) && !is.null(model$meanFunction))
    m = model$y - modelOut(model$meanFunction, model$X)
  else
    m = model$y

  ## Remove bias and apply scale.
  for (i in 1:model$d) {
    m[,i] = m[,i] - model$bias[i]
    if model$scale[i]
      m[,i] = m[,i]/model$scale[i]
  }

  return(m)
}