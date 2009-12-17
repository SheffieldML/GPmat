kernCreate <- function(x, kernType) {
  if ( is.list(x) ) {
    dim <- array()
    for ( i in 1:length(x) ) {
      dim[i] <- dim(as.matrix(x[[i]]))[2]
      if ( (dim[i] == 1) & (dim(as.matrix(x[[i]]))[1] == 1) )
          dim[i] <- x[[i]]
    }
  } else {
    dim <- dim(as.matrix(x))[2]
    if ( (dim == 1) & (dim(as.matrix(x))[1] == 1) )
      dim <- x
  }

  if ( is.list(kernType) && ("complete" %in% names(kernType)) ) {
    if ( kernType$complete == 1 ) {
      kern <- kernType
    }
    
  } else if ( is.list(kernType) ) {
    
    kern <- list(inputDimension=dim, type=kernType$type)

    if ("options" %in% names(kernType))
      kern$options <- kernType$options
    
    start <- 1    
    
    if ( kern$type == "multi" ) {
      for ( i in start:length(kernType$comp) ) {
        if ( is.list(kernType$comp) ) {
          iType <- kernType$comp[[i]]
        } else {
          iType <- kernType$comp[i]
        }
        
        if ( is.list(x) ) {
          kern$comp[[i-start+1]] <- kernCreate(x[[i-start+1]], iType)
          kern$diagBlockDim[i-start+1] <- length(x[[i-start+1]])
        } else {
          kern$comp[[i-start+1]] <- kernCreate(x, iType)
        }

        kern$comp[[i-start+1]]$index = array()
      }
        
    } else if ( kern$type %in% c("cmpnd", "tensor", "translate") )  {
      for ( i in start:length(kernType$comp) ) {
        if ( is.list(kernType$comp) ) {
          iType <- kernType$comp[[i]]
        } else {
          iType <- kernType$comp[i]
        }
        
        kern$comp[[i-start+1]] <- kernCreate(x, iType)
        kern$comp[[i-start+1]]$index = array()
      }
      
    } else if ( kern$type == "exp" ) {
      ## need double check
      if ( start == length(kernType$comp) ) {
        kern$argument <- kernCreate(x, kernType$comp[start])
      } else {
        kern$argument <- kernCreate(x, kernType$comp[start:length(kernType$comp)])
      }
    }

    kern <- kernParamInit(kern)

  } else {
    kern <- list(type=kernType, inputDimension=dim)
    if ("options" %in% names(kernType))
      kern$options <- kernType$options

    kern <- kernParamInit(kern)
  }

  kern$Kstore <- matrix()
  kern$diagK <- matrix()      

  return (kern)
  
}



kernParamInit <- function (kern) {
  
  funcName <- paste(kern$type, "KernParamInit", sep="")
  kern$transforms = list()

  func <- get(funcName, mode="function")
  kern <- func(kern)  

  return (kern)
}



kernExtractParam <- function (kern, option=1) {
  ## option=1: only return parameter values;
  ## option=2: return both parameter values and names.

  funcName <- paste(kern$type, "KernExtractParam", sep="")
  func <- get(funcName, mode="function")

  if ( option==1 ) {
    params <- func(kern)
    if ( any(is.nan(params)) )
      warning("Parameter has gone to NaN.")

    if ( "transforms" %in% names(kern) ) 
      if ( length(kern$transforms) > 0 )
        for ( i in seq(along=kern$transforms$index) ) {
          index <- kern$transforms$index[i]
          funcName <- paste(optimiDefaultConstraint(kern$transforms$type), "Transform", sep="")
          func <- get(funcName, mode="function")
          params[index] <- func(params[index], "xtoa")
        }
    
  } else {
    params <- func(kern, option)
    if ( any(is.nan(params$values)) )
      warning("Parameter has gone to NaN.")

    if ( "transforms" %in% names(kern) ) 
      if ( length(kern$transforms) > 0 )
        for ( i in seq(along=kern$transforms$index) ) {
          index <- kern$transforms$index[i]
          funcName <- paste(optimiDefaultConstraint(kern$transforms$type), "Transform", sep="")
          func <- get(funcName, mode="function")
          params$values[index] <- func(params$values[index], "xtoa")
        }             
  }

  return (params)
  
}



kernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values
  
  if ( "transforms" %in% names(kern) ) 
    if ( length(kern$transforms) > 0 )
      for ( i in seq(along=kern$transforms$index) ) {
        index <- kern$transforms$index[i]
        funcName <- paste(optimiDefaultConstraint(kern$transforms$type), "Transform", sep="")
        func <- get(funcName, mode="function")
        params[index] <- func(params[index], "atox")
      }

  funcName <- paste(kern$type, "KernExpandParam", sep="")
  func <- get(funcName, mode="function")
  kern <- func(kern, params)

  return (kern)
  
}



kernCompute <- function (kern, x, x2) {

  funcName <- paste(kern$type, "KernCompute", sep="")
  func <- get(funcName, mode="function")

  if ( nargs() < 3 ) {
    k <- func(kern, x)
  } else {
    k <- func(kern, x, x2)
  }

  return (k)
}



kernGradient <- function (kern, x, ...) {
  funcName <- paste(kern$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")

  g <- func(kern, x, ...)

  factors <- kernFactors(kern, "gradfact")
  g[factors$index] <- g[factors$index]*factors$val
  return (g)
}



kernFactors <- function (kern, factorType) {
  factors <- list(index=c(), val=c())

  if ( length(kern$transforms) > 0 ) {
    funcName <- paste(kern$type, "KernExtractParam", sep="")
    func <- get(funcName, mode="function")
    params <- func(kern)

    factors$index <- kern$transforms$index
    funcName <- paste(optimiDefaultConstraint(kern$transforms$type), "Transform", sep="")
    func <- get(funcName, mode="function")
    factors$val <- func(params[factors$index], factorType)
  }
  return (factors)
}



kernDiagCompute <- function (kern, x) {
  funcName <- paste(kern$type, "KernDiagCompute", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x)
  return (k)
}



multiKernParamInit <- function (kern) {

  kern$nParams <- 0
  kern$transforms <- list()

  if ( !("comp" %in% names(kern)) )
    kern$comp <- list()

  kern$numBlocks <- length(kern$comp)
  kern$isStationary <- TRUE

  kern$block <- list()
  for ( i in seq(along=kern$comp) ) {
    if ( !kern$comp[[i]]$isStationary )
      kern$isStationary <- FALSE

    kern$comp[[i]] <- kernParamInit(kern$comp[[i]])
    kern$nParams <- kern$nParams + kern$comp[[i]]$nParams
    kern$comp[[i]]$index <- array()

    kern$block[[i]] <- list(cross=array(), transpose=array())

    for ( j in seq(length.out=i-1) ) {
      func <- paste(kern$comp[[i]]$type, "X", kern$comp[[j]]$type, "KernCompute", sep="")
      if ( exists(func, mode="function") ) {
        kern$block[[i]]$cross[j] <- paste(kern$comp[[i]]$type, "X", kern$comp[[j]]$type, sep="")
        kern$block[[i]]$transpose[j] <- FALSE
      } else {
        func <- paste(kern$comp[[j]]$type, "X", kern$comp[[i]]$type, "KernCompute", sep="")
        if ( exists(func, mode="function") ) {
          kern$block[[i]]$cross[j] <- paste(kern$comp[[j]]$type, "X", kern$comp[[i]]$type, sep="")
          kern$block[[i]]$transpose[j] <- TRUE
        } else {
          warning(paste("No cross covariance found between", kern$comp[[i]]$type, "and", kern$comp[[j]]$type, "assuming independence."))
          kern$block[[i]]$cross[j] <- ""
          kern$block[[i]]$transpose[j] <- 0
        }
      }
    }
  }

  kern$paramGroups <- diag(1, nrow=kern$nParams, ncol=kern$nParams)

  return (kern)
}



multiKernExtractParam <- function (kern, option=1) {
  params <- cmpndKernExtractParam(kern, option)
  return (params)
}



multiKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  params <- params%*%t(kern$paramGroups)
  startVal <- 1
  endVal <- 0

  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal+kern$comp[[i]]$nParams
    kern$comp[[i]] <- kernExpandParam(kern$comp[[i]], params[startVal:endVal])
    startVal <- endVal+1
  }

  return (kern)
}



multiKernCompute <- function (kern, x, x2=x) {
  if ( is.list(x) ) {
    if ( length(x) != kern$numBlocks )
      stop ("Time information is not matched among blocks!")

    dim1 <- array(0, dim=length(x))
    dim2 <- array(0, dim=length(x))
    
    for ( i in seq(length=kern$numBlocks) ) {
      dim1[i] <- length(x[[i]])
      if ( nargs()>2 ) {
        if ( length(x) != length(x2) )
          stop ("Time information is not matched within the block!")
        dim2[i] <- length(x2[[i]])
      } else {
        dim2[i] <- dim1[i]
      }
    }

    K <- matrix(0, sum(dim1), sum(dim2))

    for ( i in seq(length=kern$numBlocks) ) {
      startOne <- sum(dim1[seq(length.out=(i-1))])+1
      endOne <- sum(dim1[seq(length.out=i)])
      startThree <- sum(dim2[seq(length.out=(i-1))])+1
      endThree <- sum(dim2[seq(length.out=i)])

      if ( nargs()<3 ) {
        K[startOne:endOne, startThree:endThree] <- multiKernComputeBlock(kern, x[[i]], i, i)
      } else {
        K[startOne:endOne, startThree:endThree] <- multiKernComputeBlock(kern, x[[i]], x2[[i]], i, i)
      }

      for ( j in seq(length.out=(i-1)) )
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          if ( nargs()<3 ) {
            K[startOne:endOne, startTwo:endTwo] <- multiKernComputeBlock(kern, x[[i]], x[[j]], i, j)
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            K[startOne:endOne, startTwo:endTwo] <- multiKernComputeBlock(kern, x[[i]], x2[[j]], i, j)
            startFour <- sum(dim1[seq(length.out=(j-1))])+1
            endFour <- sum(dim1[seq(length.out=j)])
            K[startFour:endFour, startThree:endThree] <- t(multiKernComputeBlock(kern, x2[[i]], x[[j]], j, i))
          }
        }
    }
  } else {
    # non-cell part
    dim1 = length(x)
    
    if ( nargs() > 2 ) {
      dim2 = length(x2)
    } else {
      dim2 = dim1;
    }
  
    K <- matrix(0, kern$numBlocks*dim1, kern$numBlocks*dim2)
  
    for ( i in seq(length=kern$numBlocks) ) {
      startOne <- (i-1)*dim1 + 1
      endOne <- i*dim1
      startThree <- (i-1)*dim2 + 1
      endThree <- i*dim2
      if ( nargs() < 3 ) {
        K[startOne:endOne, startThree:endThree] <- multiKernComputeBlock(kern, x, i, i)
      } else {
        K[startOne:endOne, startThree:endThree] <- multiKernComputeBlock(kern, x, x2, i, i)
      }

      for ( j in seq(length=(i-1)) ) {
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2
          if ( nargs() < 3 ) {
            K[startOne:endOne, startTwo:endTwo] <- multiKernComputeBlock(kern, x, i, j)
          } else {
            K[startOne:endOne, startTwo:endTwo] <- multiKernComputeBlock(kern, x, x2, i, j)
          }
          if ( nargs()< 3 ) {
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            startFour <- (j-1)*dim1 + 1
            endFour <- j*dim1
            K[startFour:endFour, startThree:endThree] <- t(multiKernComputeBlock(kern, x2, x, j, i))
          }
        }
      }
    }
  }

  return (K)
}



multiKernComputeBlock <- function (kern, x, x2, i, j) {
  if ( nargs() < 5 ) {
    j <- i
    i <- x2
    x2 <- array()
  }

  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernCompute", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]

    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      K <- func(arg1, x)
    } else {
      K <- func(arg1, x, x2)
    }
  } else {

    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernCompute", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernCompute", sep="")
      transpose <- !kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      arg2 <- kern$comp[[i]]
    } else {
      arg1 <- kern$comp[[i]]
      arg2 <- kern$comp[[j]]      
    }
    
    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      K <- func(arg1, arg2, x)
    } else {
      K <- func(arg1, arg2, x, x2)
    }
  }
  return (K)
}



multiKernGradient <- function (kern, x, x2, covGrad) {
  if ( is.list(x) ) {
    if ( (nargs()>3) & !is.list(x2) )
      stop("Time course information is not matched in List format.")

    arg1 <- list()
    arg2 <- list()
    dim1 <- array()
    dim2 <- array()
    for ( i in seq(length=kern$numBlocks) ) {
      dim1[i] <- length(x[[i]])
      arg1[[i]] <- x[[i]]
      if ( nargs()>3 ) {
        dim2[i] <- length(x2[[i]])
        arg2[[i]] <- x2[[i]]
      } else {
        dim2[i] <- dim1[i]
        covGrad <- x2
        arg2[[i]] <- arg1[[i]]
      }
    }

    g <- array(0, dim(kern$paramGroups)[1])
    startVal <- 1
    endVal <- 0

    for ( i in seq(length=kern$numBlocks) ) {
      endVal <- endVal + kern$comp[[i]]$nParams
      startOne <- sum(dim1[seq(length.out=(i-1))])+1
      endOne <- sum(dim1[seq(length.out=i)])
      startThree <- sum(dim2[seq(length.out=(i-1))])+1
      endThree <- sum(dim2[seq(length.out=i)])

      if ( nargs()>3 ) {
        g[startVal:endVal] <- multiKernGradientBlock(kern, arg1[[i]], arg2[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      } else {
        g[startVal:endVal] <- multiKernGradientBlock(kern, arg1[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      }

      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length.out=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          gList <- multiKernGradientBlock(kern, arg1[[i]], arg2[[j]], covGrad[startOne:endOne, startTwo:endTwo], i, j)

          g[startVal:endVal] <- g[startVal:endVal] + 2*gList$g1
          g[startVal2:endVal2] <- g[startVal2:endVal2] + 2*gList$g2
        }
        startVal2 <- endVal2 + 1
      }
      startVal <- endVal + 1
    }

  } else {
    dim1 <- length(x)
    arg1 <- x
    if ( nargs() > 3 ) {
      dim2 <- length(x2)
      arg2 <- x2
    } else {
      dim2 <- dim1
      covGrad <- x2
      arg2 <- arg1
    }

    g <- array(0, dim(kern$paramGroups)[1])
    startVal <- 1
    endVal <- 0
    for ( i in seq(length=kern$numBlocks) ) {
      endVal <- endVal + kern$comp[[i]]$nParams
      startOne <- (i-1)*dim1 + 1
      endOne <- i*dim1
      if ( nargs() > 3 ) {
        g[startVal:endVal] <- multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      } else {
        g[startVal:endVal] <- multiKernGradientBlock(kern, arg1, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      }
      
      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2

          gList <- multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, startTwo:endTwo], i, j)

          g1 <- gList$g1
          g2 <- gList$g2

          if ( nargs() > 3 ) {
            startThree <- (j-1)*dim1 + 1
            endThree <- j*dim1
            gList <- multiKernGradientBlock(kern, arg2, arg1, t(covGrad[startThree:endThree, startTwo:endTwo]), j, i)

            g3 <- gList$g1
            g4 <- gList$g2
            g[startVal:endVal] <- g[startVal:endVal] + g1 + g4
            g[startVal2:endVal2] <- g[startVal2:endVal2] + g2 + g3
          } else {
            g[startVal:endVal] <- g[startVal:endVal] + 2*g1
            g[startVal2:endVal2] <- g[startVal2:endVal2] + 2*g2           
          }
        }
        startVal2 <- endVal2 + 1
      }
      startVal <- endVal + 1
    }
  }

  g <- (g %*% kern$paramGroups)[1,]
  return (g)
}



multiKernGradientBlock <- function (kern, x, x2, covGrad, i, j) {
  if ( nargs()<6 ) {
    j <- i
    i <- covGrad
    covGrad <- x2
    x2 <- array()
  }

  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradient", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]
    factors <- kernFactors(kern$comp[[i]], "gradfact")

    func <- get(funcName, mode="function")

    if ( is.na(x2) ) {
      g <- func(arg1, x, covGrad)
    } else {
      g <- func(arg1, x, x2, covGrad)
    }
    g[factors$index] <- g[factors$index]*factors$val
    
  } else {
    if ( j<i ) {
      funcName <- paste(kern$block[[i]]$cross[j], "KernGradient", sep="")
      transpose <- kern$block[[i]]$transpose[j]
    } else {
      funcName <- paste(kern$block[[j]]$cross[i], "KernGradient", sep="")
      transpose <- kern$block[[j]]$transpose[i]
    }

    if ( transpose ) {
      arg1 <- kern$comp[[j]]
      factors1 <- kernFactors(kern$comp[[j]], "gradfact")
      arg2 <- kern$comp[[i]]
      factors2 <- kernFactors(kern$comp[[i]], "gradfact")
    } else {
      arg1 <- kern$comp[[i]]
      factors1 <- kernFactors(kern$comp[[i]], "gradfact")      
      arg2 <- kern$comp[[j]]
      factors2 <- kernFactors(kern$comp[[j]], "gradfact")
    }

    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      gList <- func(arg1, arg2, x, covGrad)
    } else {
      gList <- func(arg1, arg2, x, x2, covGrad)
    }

    g1 <- gList$g1
    g2 <- gList$g2
    
    g1[factors1$index] <- g1[factors1$index]*factors1$val
    g2[factors2$index] <- g2[factors2$index]*factors2$val

    if ( transpose ) {
      g <- g2
      g2 <- g1
      g1 <- g
    }
    g <- list(g1=g1, g2=g2)   
  
  }
  return (g)
}



multiKernDiagCompute <- function (kern, x) {
  if ( is.list(x) ) {
    dim <- 0
    for ( i in seq(along=x) )
      dim <- dim + length(x[[i]])

    k <- matrix(0, dim, 1)
    startVal <- 1
    endVal <- length(x[[1]])
    for ( i in seq(along=kern$comp) ) {

      k[startVal:endVal] <- kernDiagCompute(kern$comp[[i]], x[[i]])
      startVal <- endVal + 1
      if ( (i+1)<=length(kern$comp) )
        endVal <- endVal + length(x[[i+1]])
    }
  } else {
    k <- array(0, length(x)*kern$numBlocks)
    startVal <- 1
    endVal <- length(x)
    for ( i in seq(along=kern$comp) ) {
      k[startVal:endVal] <- kernDiagCompute(kern$comp[[i]], x)
      startVal <- endVal + 1
      endVal <- endVal + length(x)
    }
  }

  return (k)
}


multiKernCacheBlock <- function(kern, fhandle = NULL, arg = NULL, i = NULL, j = NULL, X = NULL, X2 = NULL) {

  if(!exists("cacheEnvir") || !exists("cache")) {
    cacheEnvir <- new.env("cacheEnvir")
    cache <- list()
    cache$kern$uuid <- array(dim = kern$numBlocks)
    assign("cache", cache, pos = cacheEnvir)
    K <- array()
    return (K)
  }

  cache <- get("cache", envir = as.environment(cacheEnvir)) 
  myCache <- cache$kern$uuid[i,j]
  key <- rbind(X, X2)

  for (k in 1:length(myCache)) {
    if (dim(key)[2] == length(myCache[k][1]) && all(key == myCache[k][1])) {
      K = myCache[k][2]
      return (K)
    }
  }

  # No match if we get all the way here
  if (is.null(X2)) {
    K = fhandle(arg, X)
  }
  else {
    K = fhandle(arg, X, X2)
  }
  cache$kern$uuid[i, j][end+1] <- c(key, K)
}



multiKernFixBlocks <- function(kern, blocks = 1) {

  length <- 30
  op <- options(digits.secs=15)
  time <- Sys.time()
  uuid <- ""
  for (i in 1:length) {
    char <- substr(time, i, i)
    if (is.na(charmatch(char, ".")) && is.na(charmatch(char, ":")) && is.na(charmatch(char, "-")) && is.na(charmatch(char, " "))) {
      uuid <- paste(uuid, char, sep="")
    }
  }
  uuid <- paste('a', uuid, sep="")

  kern$fixedBlocks <- array(0, dim = c(1, kern$numBlocks))
  kern$fixedBlocks[blocks] <- 1
  kern$uuid <- uuid
  multiKernCacheBlock(kern)

  return(kern)
}
