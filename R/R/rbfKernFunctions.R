rbfKernParamInit <- function (kern) {
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$nParams <- 2
  kern$paramNames <- c("inverseWidth", "variance")
  
  kern$isStationary <- TRUE

  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && kern$options$isNormalised)
    kern$isNormalised <- TRUE
  else
    kern$isNormalised <- FALSE

  if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options)) {
    kern$transforms <- list(list(index=1, type="bounded"),
                            list(index=2, type="positive"))
    kern$transformArgs <- list()
    kern$transformArgs[[1]] <- kern$options$inverseWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  }
  else {
    kern$transforms <- list(list(index=c(1,2), type="positive"))
  }

  return (kern)
}



# untransformed.values is ignored
rbfKernExtractParam <- function (kern, only.values=TRUE,
                                 untransformed.values=TRUE) {
  params <- c(kern$inverseWidth, kern$variance)

  if ( !only.values )
    names(params) <- c("inverseWidth", "variance")

  return (params)
}



rbfKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]
  kern$variance <- params[2]

  return (kern)
}


rbfKernDisplay <- function (kern, spaceNum=0) {
  spacing <- matrix("", spaceNum+1)
  cat(spacing)
  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    cat("Normalised version of the kernel.\n")
  else
    cat("Unnormalised version of the kernel.\n")
  cat(spacing)
  cat("RBF inverse width: ", kern$inverseWidth, " (length scale ", 1/sqrt(kern$inverseWidth), ")\n", sep="")
  cat(spacing)
  cat("RBF variance: ", kern$variance, "\n", sep="")
}

rbfKernCompute <- function (kern, x, x2=NULL) {
  if ( nargs() < 3 ) {
    n2 <- .dist2(x,x)
  } else {
    n2 <- .dist2(x,x2)
  }

  wi2 <- 0.5*kern$inverseWidth
  k <- kern$variance*exp(-n2*wi2)

  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))
  
  return (k)
}



rbfKernGradient <- function (kern, x, x2, covGrad) {
  if ( nargs()==3 ) {
    k <- rbfKernCompute(kern, x)
    dist2xx <- .dist2(x, x)
    covGrad <- x2
  } else if ( nargs()==4 ) {
    k <- rbfKernCompute(kern, x, x2)
    dist2xx <- .dist2(x, x2)
  }

  g <- array()
  if ("isNormalised" %in% names(kern) && kern$isNormalised) {
    g[1] <- -0.5*sum(covGrad*k*dist2xx) +
      0.5 * sum(covGrad*k)/kern$inverseWidth
  }
  else {
    g[1] <- -0.5*sum(covGrad*k*dist2xx)
  }
  g[2] <- sum(covGrad*k)/kern$variance

  if ( any(is.nan(g)) )
    warning("g is NaN.")

  return (g)
}



rbfKernDiagCompute <- function (kern, x) {
  k <- matrix(kern$variance, dim(as.array(x))[1], 1)

  if ("isNormalised" %in% names(kern) && kern$isNormalised)
    k <- k * sqrt(kern$inverseWidth/(2*pi))

  return (k)
}
