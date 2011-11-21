modelOutputGrad <- function(model, X, dim) {

  if (nargs() > 2) {
    fhandle = get(paste(model$type, "OutputGrad", sep=""), mode="function")
    g = fhandle(model, X, dim)
  } else {
    fhandle = get(paste(model$type, "OutputGrad", sep=""), mode="function")
    gtemp = fhandle(model, X)
    if ("paramGroups" %in% names(model)) {
      g = matrix(0, dim(X)[1], dim(model$paramGroups)[2], dim(gtemp)[3])
      for (i in 1:dim(gtemp)[3])
	g = gtemp[, , i]%*%model$paramGroups
    } else
      g = gtemp
  }

  return (g)
}
