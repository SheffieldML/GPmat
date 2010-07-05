mlpKernParamInit <- function (kern) {
  kern$weightVariance <- 10
  kern$biasVariance <- 10
  kern$variance <- 1
  kern$nParams <- 3

  kern$transforms <- list(list(index=c(1,2,3), type="positive"))

  kern$isStationary <- FALSE

  return (kern)
}



# untransformed.values is ignored
mlpKernExtractParam <- function (kern, only.values=TRUE,
                                 untransformed.values=TRUE) {
  params <- c(kern$weightVariance, kern$biasVariance, kern$variance)

  if ( !only.values )
    names(params) <- c("weightVariance", "biasVariance", "variance")

  if ( any(is.na(params)) )
    warning("params are infinite")

  return (params)
}



mlpKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  if ( any(is.na(params)) )
    warning("params are infinite")

  kern$weightVariance <- params[1]
  kern$biasVariance <- params[2]
  kern$variance <- params[3]

  return (kern)
}



mlpKernCompute <- function (kern, x, x2=NULL, option=1) {
 
  if ( (nargs()<3) | (nargs()==4 & is.null(x2)) ) {
    innerProd <- x %*% t(x)
    numer <- innerProd*kern$weightVariance + kern$biasVariance
    vec <- diag(numer) + 1
    denom <- sqrt(vec %*% t(vec))
    arg <- numer/denom
    k <- kern$variance*asin(arg)
  } else {
    innerProd <- x %*% t(x2)
    numer <- innerProd*kern$weightVariance + kern$biasVariance
    vec1 <- apply(as.matrix(x)*as.matrix(x), 1, sum) * kern$weightVariance + kern$biasVariance + 1
    vec2 <- apply(as.matrix(x2)*as.matrix(x2), 1, sum) * kern$weightVariance + kern$biasVariance + 1
    denom <- sqrt(vec1 %*% t(vec2))
    arg <- numer/denom
    k <- kern$variance*asin(arg)
  }

  if ( option == 2 )
    k <- list(k=k, innerProd=innerProd, arg=arg, denom=denom, numer=numer)

  return (k)
}



mlpKernGradient <- function (kern, x, x2, covGrad) {

  if ( nargs()==3 ) {
    mlpK <- mlpKernCompute(kern, x, NULL, 2)
    k <- mlpK$k
    innerProd <- mlpK$innerProd
    arg <- mlpK$arg
    denom <- mlpK$denom
    numer <- mlpK$numer
    covGrad <- x2
  } else if ( nargs()==4 ) {
    mlpK <- mlpKernCompute(kern, x, x2, 2)
    k <- mlpK$k
    innerProd <- mlpK$innerProd
    arg <- mlpK$arg
    denom <- mlpK$denom
    numer <- mlpK$numer
  }

  denom3 <- denom*denom*denom
  base <- kern$variance/sqrt(1-arg*arg)
  baseCovGrad <- base*covGrad

  g <- array(0,3)
  if ( nargs()==3 ) {
    vec <- diag(innerProd)
    g[1] <- sum((innerProd/denom-0.5*numer/denom3*((kern$weightVariance*vec+kern$biasVariance+1)%*%t(vec)+vec%*%t(kern$weightVariance*vec+kern$biasVariance+1)))*baseCovGrad)
    g[2] <- sum((1/denom-0.5*numer/denom3*(matrix(rep(vec, times=length(vec)), length(vec), length(vec))*kern$weightVariance+2*kern$biasVariance+2+matrix(rep(vec, each=length(vec)), length(vec), length(vec))*kern$weightVariance))*baseCovGrad)
  } else {
    vec1 <- apply(x*x, 1, sum)
    vec2 <- apply(x2*x2, 1, sum)
    g[1] <- sum((innerProd/denom-0.5*numer/denom3*((kern$weightVariance*vec1+kern$biasVariance+1)*t(vec2)+vec1*t(kern$weightVariance*vec2+kern$biasVariance+1)))*baseCovGrad)
    g[2] <- sum((1/denom-0.5*numer/denom3*(t(matrix(rep(vec1, each=length(vec2)), length(vec2), length(vec1)))*kern$weightVariance+2*kern$biasVariance+2+matrix(rep(vec2, each=length(vec1)), length(vec1), length(vec2))*kern$weightVariance))*baseCovGrad)
  }

  g[3] <- sum(k/kern$variance*covGrad)

  if ( any(is.na(g)) )
    warning("g is NA.")

  return (g)
}



mlpKernDiagGradX <- function (kern, X) {
  gX <- array(0, dim=dim(as.matrix(X)))
  for ( i in seq(length=dim(X)[1]) ) {
    gX[i,] <- .mlpKernDiagGradXpoint(kern, as.matrix(X)[i,])
  }

  return (gX)
}



.mlpKernDiagGradXpoint <- function (kern, x) {
  x <- as.matrix(x)
  innerProd <- x %*% t(x)
  numer <- innerProd*kern$weightVariance + kern$biasVariance
  demon <- numer + 1
  arg <- numer / demon
  gX <- array(0, dim = dim(x))
  for ( j in seq(length=dim(x)[2]) ) {
    gX[,j] <- 1/demon - numer/demon^2
    gX[,j] <- ((2*x[,j]*kern$weightVariance*kern$variance) %*% gX[,j] )/sqrt(1-arg*arg)
  }

  return (gX)
}



mlpKernGradX <- function (kern, X, X2) {
  gX <- array(0, dim=c(dim(as.matrix(X2)), dim(X)[1]))
  for ( i in seq(length=dim(X)[1]) ) {
    gX[,,i] <- .mlpKernGradXpoint(kern, as.matrix(X)[i,], X2)
  }

  return (gX)
}



.mlpKernGradXpoint <- function (kern, x, X2) {
  x <- as.matrix(x)
  X2 <- as.matrix(X2)
  innerProd <- X2 %*% t(x)
  numer <- innerProd*kern$weightVariance + kern$biasVariance
  vec1 <- apply(x*x, 1, sum)*kern$weightVariance + kern$biasVariance + 1
  vec2 <- apply(X2*X2, 1, sum)*kern$weightVariance + kern$biasVariance + 1
  demon <- sqrt(vec2 %*% t(vec1))
  arg <- numer / demon
  gX <- array(0, dim = dim(X2))
  for ( j in seq(length=dim(X2)[2]) ) {
    gX[,j] <- X2[,j]/demon - vec2*x[,j]*numer/demon^3
    gX[,j] <- kern$weightVariance*kern$variance*gX[,j]/sqrt(1-arg*arg)
  }

  return (gX)
}
