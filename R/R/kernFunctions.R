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

  if ( is.list(kernType) && kernType$type == "parametric" ) {
    kernOptions <- kernType$options
    kernType <- kernType$realType
  }
  else
    kernOptions <- NULL
  
  if ( is.list(kernType) && ("complete" %in% names(kernType)) ) {
    if ( kernType$complete == 1 ) {
      kern <- kernType
    }
    
  } else if ( is.list(kernType) ) {
    
    kern <- list(inputDimension=dim, type=kernType$type)

    if (!is.null(kernOptions))
      kern$options <- kernOptions
    
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

    if (!is.null(kernOptions))
      kern$options <- kernOptions

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



kernExtractParam <- function (kern, only.values=TRUE) {
  funcName <- paste(kern$type, "KernExtractParam", sep="")
  func <- get(funcName, mode="function")

  params <- func(kern, only.values=only.values)

  if ( any(is.nan(params)) )
    warning("Parameter has gone to NaN.")

  if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0) )
    for ( i in seq(along=kern$transforms) ) {
      index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "xtoa", kern$transformArgs[[i]])
      else
        params[index] <- func(params[index], "xtoa")
    }

  return (params)
}



kernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values
  
  if ( "transforms" %in% names(kern) && (length(kern$transforms) > 0) )
    for ( i in seq(along=kern$transforms) ) {
      index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "atox", kern$transformArgs[[i]])
      else
        params[index] <- func(params[index], "atox")
    }

  funcName <- paste(kern$type, "KernExpandParam", sep="")
  func <- get(funcName, mode="function")
  kern <- func(kern, params)

  return (kern)
  
}

kernDisplay <- function (kern, ...) {

  funcName <- paste(kern$type, "KernDisplay", sep="")
  if(exists(funcName, mode="function")) {
    func <- get(funcName, mode="function")
    return (func(kern, ...))
  }

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

  factors <- .kernFactors(kern, "gradfact")
  for (i in seq(along=factors))
    g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val

  return (g)
}



.kernFactors <- function (kern, factorType) {
  factors <- list()

  if ( length(kern$transforms) > 0 ) {
    funcName <- paste(kern$type, "KernExtractParam", sep="")
    func <- get(funcName, mode="function")
    params <- func(kern)

    for (i in seq(along=kern$transforms)) {
      factors[[i]] <- list()
      factors[[i]]$index <- kern$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(kern$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType, kern$transformArgs[[i]])
      else
        factors[[i]]$val <- func(params[factors[[i]]$index], factorType)
    }
  }
  return (factors)
}



kernDiagCompute <- function (kern, x) {
  funcName <- paste(kern$type, "KernDiagCompute", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x)
  return (k)
}



kernDiagGradX <- function (kern, x) {
  funcName <- paste(kern$type, "KernDiagGradX", sep="")
  func <- get(funcName, mode="function")
  k <- func(kern, x)
  return (k)
}



kernGradX <- function (kern, x) {
  funcName <- paste(kern$type, "KernGradX", sep="")
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

  kern$fixBlocks <- rep(FALSE, kern$numBlocks)
  if ("options" %in% names(kern) && "fixedBlocks" %in% names(kern$options)
      && kern$options$fixedBlocks) {
    kern$fixBlocks[kern$options$fixedBlocks] <- TRUE
    kern$cache <- new.env(parent=emptyenv())
    assign("cache", list(list(list())), envir=kern$cache)
  }

  return (kern)
}



multiKernExtractParam <- function (kern, only.values=TRUE) {
  return (cmpndKernExtractParam(kern, only.values=only.values))
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

multiKernDisplay <- function (kern, spaceNum=0) {
  spacing = matrix("", spaceNum+1)
  cat(spacing);
  cat("Multiple output block kernel:\n")
  for(i in seq(along=kern$comp)) {
    cat(spacing)
    cat("Block ", i, "\n", sep="")
    kernDisplay(kern$comp[[i]], spaceNum)
  }
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
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x[[i]])
      } else {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x[[i]], x2[[i]])
      }

      for ( j in seq(length.out=(i-1)) )
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          if ( nargs()<3 ) {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x[[i]], x[[j]])
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x[[i]], x2[[j]])
            startFour <- sum(dim1[seq(length.out=(j-1))])+1
            endFour <- sum(dim1[seq(length.out=j)])
            K[startFour:endFour, startThree:endThree] <- t(.multiKernComputeBlock(kern, j, i, x2[[i]], x[[j]]))
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
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x)
      } else {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x, x2)
      }

      for ( j in seq(length=(i-1)) ) {
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2
          if ( nargs() < 3 ) {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x)
          } else {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x, x2)
          }
          if ( nargs()< 3 ) {
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            startFour <- (j-1)*dim1 + 1
            endFour <- j*dim1
            K[startFour:endFour, startThree:endThree] <- t(.multiKernComputeBlock(kern, j, i, x2, x))
          }
        }
      }
    }
  }

  return (K)
}



.multiKernComputeBlock <- function (kern, i, j, x1, x2=NULL) {
  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernCompute", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]

    func <- get(funcName, mode="function")
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, x1)
      else
        K <- func(arg1, x1, x2)
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
    if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
      K <- .multiKernCacheBlock(kern, func, i, j, arg1=arg1, arg2=arg2, x1=x1, x2=x2)
    }
    else {
      if (is.null(x2))
        K <- func(arg1, arg2, x1)
      else
        K <- func(arg1, arg2, x1, x2)
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
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1[[i]], arg2[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      } else {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      }

      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length.out=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          gList <- .multiKernGradientBlock(kern, arg1[[i]], arg2[[j]], covGrad[startOne:endOne, startTwo:endTwo], i, j)

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
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      } else {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      }
      
      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2

          gList <- .multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, startTwo:endTwo], i, j)

          g1 <- gList$g1
          g2 <- gList$g2

          if ( nargs() > 3 ) {
            startThree <- (j-1)*dim1 + 1
            endThree <- j*dim1
            gList <- .multiKernGradientBlock(kern, arg2, arg1, t(covGrad[startThree:endThree, startTwo:endTwo]), j, i)

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



.multiKernGradientBlock <- function (kern, x, x2, covGrad, i, j) {
  if ( nargs()<6 ) {
    j <- i
    i <- covGrad
    covGrad <- x2
    x2 <- array()
  }

  if (kern$fixBlocks[[i]] && kern$fixBlocks[[j]]) {
    if (i==j)
      return (0)
    else
      return (list(g1=0, g2=0))
  }

  if ( i==j ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradient", sep="")
    transpose <- 0
    arg1 <- kern$comp[[i]]
    factors <- .kernFactors(kern$comp[[i]], "gradfact")

    func <- get(funcName, mode="function")

    if ( is.na(x2) ) {
      g <- func(arg1, x, covGrad)
    } else {
      g <- func(arg1, x, x2, covGrad)
    }
    for (i in seq(along=factors))
      g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val
    
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
      factors1 <- .kernFactors(kern$comp[[j]], "gradfact")
      arg2 <- kern$comp[[i]]
      factors2 <- .kernFactors(kern$comp[[i]], "gradfact")
    } else {
      arg1 <- kern$comp[[i]]
      factors1 <- .kernFactors(kern$comp[[i]], "gradfact")      
      arg2 <- kern$comp[[j]]
      factors2 <- .kernFactors(kern$comp[[j]], "gradfact")
    }

    func <- get(funcName, mode="function")
    if ( any(is.na(x2)) ) {
      gList <- func(arg1, arg2, x, covGrad)
    } else {
      gList <- func(arg1, arg2, x, x2, covGrad)
    }

    g1 <- gList$g1
    g2 <- gList$g2
    
    for (i in seq(along=factors1))
      g1[factors1[[i]]$index] <- g1[factors1[[i]]$index]*factors1[[i]]$val

    for (i in seq(along=factors2))
      g2[factors2[[i]]$index] <- g2[factors2[[i]]$index]*factors2[[i]]$val

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


.multiKernCacheBlock <- function(kern, fhandle, i, j, x1, x2=NULL, arg1, arg2=NULL) {

  global_cache <- get("cache", envir = kern$cache)
  if ((length(global_cache) >= i) && (length(global_cache[[i]]) >= j))
    cache <- global_cache[[i]][[j]]
  else
    cache <- list()
  key <- c(x1, x2)

  for (k in seq(along=cache)) {
    if (length(key) == length(cache[[k]]$key) && all(key == cache[[k]]$key)) {
      #cat("multiKernCacheBlock: cache hit ", i, j, x1, x2, "\n")
      return (cache[[k]]$value)
    }
  }

  #cat("multiKernCacheBlock: cache miss ", i, j, x1, x2, "\n")
  # No match if we get all the way here
  if (is.null(arg2)) {
    if (is.null(x2))
      K <- fhandle(arg1, x1)
    else
      K <- fhandle(arg1, x1, x2)
  }
  else {
    if (is.null(x2))
      K <- fhandle(arg1, arg2, x1)
    else
      K <- fhandle(arg1, arg2, x1, x2)
  }
  cache <- append(cache, list(list(key=key, value=K)))
  if (length(global_cache) < i)
    global_cache[[i]] <- list()
  global_cache[[i]][[j]] <- cache
  assign("cache", global_cache, envir=kern$cache)
  return(K)
}
