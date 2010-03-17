whiteKernParamInit <- function (kern) {
  
  kern$variance <- exp(-2)
  kern$nParams <- 1
  kern$paramNames <- c("variance")
  
  kern$transforms <- list(list(index=c(1), type="positive"))

  kern$isStationary <- TRUE

  return (kern)
}


whiteXwhiteKernCompute <- function (whiteKern1, whiteKern2, x1, x2) {
  if ( nargs()<4 )
    x2=x1

  k <- matrix(0, nrow=length(x1), ncol=length(x2))
  return (k)  
}


whiteKernDisplay <- function (kern, spaceNum=0) {
  spacing <- matrix("", spaceNum+1)
  cat(spacing)
  cat("White Noise Variance: ", kern$variance, "\n", sep="")
}


whiteKernCompute <- function (kern, x, x2) {
  if ( nargs()<3 ) {
    k <- kern$variance*diag(1, nrow=length(x), ncol=length(x))
  } else {
    k <- matrix(0, nrow=length(x), ncol=length(x2))
  }
  return (k)
}



whiteKernExtractParam <- function (kern, only.values=TRUE) {
  params <- c(kern$variance)

  if ( !only.values ) {
    names(params) <- c("variance")
  }

  return (params)
}



whiteKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$variance <- params[1]

  return (kern)
}



whiteKernGradient <- function (kern, x, x2, covGrad) {
  if ( nargs()==3 ) {
    covGrad <- x2
    g <- sum(diag(covGrad))
  } else {
    g <- 0
  }  

  return (g)
}



whiteXwhiteKernGradient <- function (whiteKern1, whiteKern2, x1, x2, covGrad) {
  if ( nargs()<5 ) {
    covGrad <- x2
    x2 <- x1
  }
  
  if ( dim(as.matrix(x1))[2]>1 | dim(as.matrix(x2))[2]>1 )
    stop("Input can only have one column.")

  g <- list(g1=0, g2=0)
  return (g)
}



whiteKernDiagCompute <- function (kern, x) {
  k <- matrix(kern$variance, length(x), 1)
  return (k)
}
