cmpndKernParamInit <- function (kern) {
  
  kern$nParams <- 0
  kern$transforms <- list()

  if ( !any(grep("comp", names(kern))) )
    kern$comp <- list()

  for ( i in seq(along=kern$comp) ) {

    kern$comp[[i]] <- kernParamInit(kern$comp[[i]])
    kern$nParams <- kern$nParams + kern$comp[[i]]$nParams
    kern$comp[[i]]$index <- array()

    if ( any(grep("numBlocks", names(kern$comp[[i]]))) ) {
      if ( i==1 ) {
        kern$numBlocks <- kern$comp[[i]]$numBlocks
      } else {
        if ( (!any(grep("numBlocks", names(kern)))) | (kern$numBlocks!=kern$comp[[i]]$numBlocks) ) {
          stop("Compound of multi kernels with different numbers of blocks.")
        }
      }
    } else {
      if ( any(grep("numBlocks", names(kern))) )
        stop("Attempt to combine multi-kernel with non multi-kernel.")
    }
  }

  kern$paramGroups <- diag(1, nrow=kern$nParams, ncol=kern$nParams)

  kern$whiteVariance <- 0
  kern$isStationary <- TRUE

  for ( i in seq(along=kern$comp) ) {
    if ( !kern$comp[[i]]$isStationary )
      kern$isStationary <- FALSE

    if ( kern$comp[[i]]$type == "white" ) {
      kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$variance
    } else {
      if ( any(grep("whiteVariance", names(kern$comp[[i]]))) ) {
        kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$whiteVariance
      }
    }
  }

  return (kern)
  
}



cmpndKernExtractParam <- function (kern, option=1) {

  startVal <- 1
  endVal <- 0
  
  if ( option==1 ) {
    params <- c()

    for ( i in seq(along=kern$comp) ) 
      params <- c(params, kernExtractParam(kern$comp[[i]]))

  } else {
    storedTypes <- c()
    params <- list(values=c(), names=c())
    origNames <- c()
    for ( i in seq(along=kern$comp) ) {
      paramsList <- kernExtractParam(kern$comp[[i]], option)
      params$values <- c(params$values, paramsList$values)
      kernName <- paste(kern$comp[[i]]$type, length(grep(kern$comp[[i]]$type, storedTypes))+1, sep="")
      paramName <- paste(kernName, paramsList$names, sep="_")
      origNames <- c(origNames, paramName)
      storedTypes <- c(storedTypes, kern$comp[[i]]$type)
    }
  }

  if ( any(grep("paramGroups", names(kern))) ) {
    paramGroups <- kern$paramGroups
    for ( i in seq(length.out=dim(paramGroups)[2]) ) {
      ind <- grep(1, paramGroups[,i])
      if ( is.list(params) ) {
        params$names[i] <- origNames[ind[1]]
        for ( j in seq(2, length.out=length(ind)-1) )
          params$names[i] <- paste(params$names[i], origNames[ind[j]],sep="/")
      }
   
      paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
    }
  }

  if ( is.list(params) ) {
    params$values <- params$values%*%paramGroups
  } else {
    params <- params%*%paramGroups
  }    

  return (params)

}



cmpndKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values
  params <- params %*% t(kern$paramGroups)
  startVal <- 1
  endVal <- 0
  kern$whiteVariance <- 0
  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal+kern$comp[[i]]$nParams
    kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[startVal:endVal])
    startVal <- endVal+1
    if ( any(grep("white", kern$comp[[i]]$type)) ) {
      kern$whiteVariance <- kern$whiteVairance+kern$comp[[i]]$variance
    } else if ( any(grep("whiteVariance", names(kern$comp[[i]]))) ) {
      kern$whiteVariance <- kern$whiteVariance+kern$comp[[i]]$whiteVariance
    }
  }      

  return (kern)
}



cmpndKernCompute <- function (kern, x, x2) {
  if ( nargs()>2 ) {
    i <- 1
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index])
    } else {
      k <- kernCompute(kern$comp[[i]], x, x2)
    }
    for ( i in seq(2, length.out=(length(kern$comp)-1)) )
      if ( !is.na(kern$comp[[i]]$index) ) {
        k <- k+kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index])
      } else {
        k <- k+kernCompute(kern$comp[[i]], x, x2)
      }
  } else {
    i <- 1
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
    } else {
      k <- kernCompute(kern$comp[[i]], x)
    }
    for ( i in seq(2, length.out=(length(kern$comp)-1)) )
      if ( !is.na(kern$comp[[i]]$index) ) {
        k <- k+kernCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
      } else {
        k <- k+kernCompute(kern$comp[[i]], x)
      }
  }
  return (k)  
}



cmpndKernGradient <- function (kern, x, x2, covGrad) {

  if ( nargs()<4 ) 
    covGrad <- x2
  
  g <- array(0, dim(kern$paramGroups)[1])
  startVal <- 1
  endVal <- 0

  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal + kern$comp[[i]]$nParams
    if ( !is.na(kern$comp[[i]]$index) ) {
      if ( nargs() < 4 ) {
        g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
      } else {
        g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
      }
    } else {
      if ( nargs() < 4 ) {
        g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, covGrad)
      } else {
        g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, x2, covGrad)
      }
    }
    startVal <- endVal + 1       
  }

  g <- g %*% kern$paramGroups    

  return (g)
}



cmpndKernDiagCompute <- function (kern, x) {
  i <- 1
  if ( !is.na(kern$comp[[i]]$index) ) {
    k <- kernDiagCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
  } else {
    k <- kernDiagCompute(kern$comp[[i]], x)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) )
    if ( !is.na(kern$comp[[i]]$index) ) {
      k <- k + kernDiagCompute(kern$comp[[i]], x[,kern$comp[[i]]$index])
    } else {
      k <- k + kernDiagCompute(kern$comp[[i]], x)
    }
       
  return (k)
}



cmpndKernGradX <- function (kern, X, X2) {
  i <- 1
  funcName <- paste(kern$comp[[i]]$type, "KernGradX", sep="")
  func <- get(funcName, mode="function")

  if ( !is.na(kern$comp[[i]]$index) ) {
    gX <- array(0, dim=c(length(X2), length(X2), length(X)))
    gX[,kern$comp[[i]]$index,] <- func(kern$comp[[i]], X[,kern$comp[[i]]$index], X2[,kern$comp[[i]]$index])
  } else {
    gX <- func(kern$comp[[i]], X, X2)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradX", sep="")
    func <- get(funcName, mode="function")
    if ( !is.na(kern$comp[[i]]$index) ) {
      gX[,kern$comp[[i]]$index,] <- gX[,kern$comp[[i]]$index,] +  func(kern$comp[[i]], X[,kern$comp[[i]]$index], X2[,kern$comp[[i]]$index])
    } else {
      gX <- gX + func(kern$comp[[i]], X, X2)
    }
  }
    
  return (gX)
}



cmpndKernDiagGradX <- function (kern, X) {
  X <- as.matrix(X)
  i <- 1
  funcName <- paste(kern$comp[[i]]$type, "KernDiagGradX", sep="")
  func <- get(funcName, mode="function")

  if ( !is.na(kern$comp[[i]]$index) ) {
    gX <- array(0, dim=dim(X))
    gX[,kern$comp[[i]]$index,] <- kernDiagGradX(kern$comp[[i]], X[,kern$comp[[i]]$index])
  } else {
    gX <- func(kern$comp[[i]], X)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) ) {
    if ( !is.na(kern$comp[[i]]$index) ) {
      gX[,kern$comp[[i]]$index] <- gX[,kern$comp[[i]]$index] +  func(kern$comp[[i]], X[,kern$comp[[i]]$index])
    } else {
      gX <- gX + func(kern$comp[[i]], X)
    }
  }
    
  return (gX)
}
