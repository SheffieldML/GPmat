translateKernDiagCompute <- function (kern, x) {
  x <- x - array(kern$centre, length(x))
  k <- cmpndKernDiagCompute(kern, x)
       
  return (k)
}



translateKernExtractParam <- function (kern, only.values=TRUE) {
  kern$nParams <- kern$nParams - kern$inputDimension
  
  params <- cmpndKernExtractParam(kern, only.values=only.values)
  centre <- c(kern$centre)
  if ( !only.values ) {
    for ( i in seq(length=kern$inputDimension) ) {
      names(centre)[i] <- paste("Centre", i, sep="")
    }
  }
  params <- c(params, centre)

  return (params)
}



translateKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  endVal <- length(params) - kern$inputDimension
  kern <- cmpndKernExpandParam(kern, params[1:endVal])
  startVal <- endVal+1
  endVal <- length(params)
  kern$centre <- params[startVal:endVal]

  return (kern)
}



translateKernCompute <- function (kern, x, x2) {
  if ( nargs()>2 ) {
    x <- x - array(kern$centre, length(x))
    x2 <- x2 - array(kern$centre, length(x2))
    k <- cmpndKernCompute(kern, x, x2)
  } else {
    x <- x - array(kern$centre, length(x))
    k <- cmpndKernCompute(kern, x)
  }
  return (k)  
}



translateKernParamInit <- function (kern) {
  kern <- cmpndKernParamInit(kern)

  kern$nParams <- kern$nParams + kern$inputDimension
  kern$centre <- array(0, kern$inputDimension)

  return (kern)
}



translateKernGradient <- function (kern, x, x2, covGrad) {

  if ( nargs()<4 ) {
    covGrad <- x2
    x <- x - array(kern$centre, length(x))
    g <- cmpndKernGradient(kern, x, covGrad)

    gKX <- cmpndKernGradX(kern, x, x)
    gKX <- gKX*2
    dgKX <- cmpndKernDiagGradX(kern, x)
    for ( i in seq(length=dim(gKX)[1]) )
      gKX[i,,i] <- dgKX[i,]
   
  } else {
    x <- x - array(kern$centre, length(x))
    x2 <- x2 - array(kern$centre, length(x2))
    g <- cmpndKernGradient(kern, x, x2, covGrad)

    gKX12 <- cmpndKernGradX(kern, x, x2)
    gKX21 <- cmpndKernGradX(kern, x2, x)
  }

  gcentre <- array(0, kern$inputDimension)
  if ( nargs()<4 ) {
    for ( i in seq(length=length(x)) ) {
      for ( j in seq(length=kern$inputDimension) ) 
        gcentre[j] <- gcentre[j] - t(gKX[,j,i]) %*% covGrad[,i]
    }
  } else {
    for ( i in seq(length=length(x)) ) {
      for ( j in seq(length=kern$inputDimension) ) 
        gcentre[j] <- gcentre[j] - sum(gKX12[,j,i]*t(covGrad[i,]))
    }
    for ( i in seq(length=length(x2)) ) {
      for ( j in seq(length=kern$inputDimension) ) 
        gcentre[j] <- gcentre[j] - sum(gKX21[,j,i]*t(covGrad[,i]))
    }
  }

  g <- c(g, gcentre)

  return (g)
}
