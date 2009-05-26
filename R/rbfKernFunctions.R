rbfKernParamInit <- function (kern) {
  kern$inverseWidth <- 1
  kern$variance <- 1
  kern$nParams <- 2

  kern$transforms <- list(index=c(1,2), type="positive")

  kern$isStationary=TRUE

  return (kern)
}



rbfKernExtractParam <- function (kern, option=1) {

  if ( option == 1 ) {
    params <- c(kern$inverseWidth, kern$variance)
    
  } else {
    params <- list(values=c(kern$inverseWidth, kern$variance), names=c("inverseWidth", "variance"))    
  }

  return (params)
}



rbfKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$inverseWidth <- params[1]
  kern$variance <- params[2]

  return (kern)
}



rbfKernCompute <- function (kern, x, x2=NULL) {
  if ( nargs() < 3 ) {
    n2 <- dist2(x,x)
  } else {
    n2 <- dist2(x,x2)
  }

  wi2 <- 0.5*kern$inverseWidth
  k <- kern$variance*exp(-n2*wi2)

  return (k)
}



rbfKernGradient <- function (kern, x, x2, covGrad) {
  if ( nargs()==3 ) {
    k <- rbfKernCompute(kern, x)
    dist2xx <- dist2(x, x)
    covGrad <- x2
  } else if ( nargs()==4 ) {
    k <- rbfKernCompute(kern, x, x2)
    dist2xx <- dist2(x, x2)
  }

  g <- array()
  g[1] <- -0.5*sum(covGrad*k*dist2xx)
  g[2] <- sum(covGrad*k)/kern$variance

  if ( any(is.nan(g)) )
    warning("g is NaN.")

  return (g)
}



rbfKernDiagCompute <- function (kern, x) {
  k <- matrix(kern$variance, length(x), 1)
  return (k)
}
