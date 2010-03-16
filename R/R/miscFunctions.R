expTransform <- function (x, transform="atox") {

  eps <- 2.2204e-16
  maxVal <- 4.3112e+15

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( "atox" == transform ) {
    for ( ind in 1:length(as.array(x)) ) {
      if ( x[ind] > thre ) y[ind] <- maxVal else
      if ( x[ind] < -thre ) y[ind]<- eps else
      y[ind] <- exp(x[ind])
    }
  } else if ( "xtoa" == transform ) {
    for ( ind in 1:length(as.array(x)) ) {
      y[ind] <- .complexLog(x[ind])
    }
  } else if ( "gradfact" == transform )
    y <- x

  return (y)
}


sigmoidTransform <- function (x, transform="atox") {

  eps <- 2.2204e-16

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( "atox" == transform ) {
    for ( ind in 1:length(as.array(x)) ) {
      if ( x[ind] > thre )
        y[ind] <- 1-eps
      else if ( x[ind] < -thre )
        y[ind]<- eps
      else
        y[ind] <- 1/(1+exp(-x[ind]))
    }
  } else if ( "xtoa" == transform ) {
    for ( ind in 1:length(as.array(x)) ) {
      y[ind] <- .complexLog(x[ind]/(1-x[ind]))
    }
  } else if ( "gradfact" == transform )
    y <- x*(1-x)

  return (y)
}


boundedTransform <- function (x, transform="atox", bounds) {

  eps <- 2.2204e-16

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( "atox" == transform ) {
    for ( ind in 1:length(as.array(x)) ) {
      if ( x[ind] > thre )
        y[ind] <- 1-eps
      else if ( x[ind] < -thre )
        y[ind] <- eps
      else
        y[ind] <- 1/(1+exp(-x[ind]))
    }
    y <- (bounds[2] - bounds[1])*y + bounds[1]
  } else if ( "xtoa" == transform ) {
    x <- (x - bounds[1]) / (bounds[2] - bounds[1])
    for ( ind in 1:length(as.array(x)) ) {
      y[ind] <- .complexLog(x[ind]/(1-x[ind]))
    }
  } else if ( "gradfact" == transform ) {
    y <- (x-bounds[1])*(1-(x-bounds[1])/(bounds[2] - bounds[1]))
  }

  return (y)
}


modelTieParam <- function (model, paramsList) {
  columnToDel <- c()

  if ( !is.list(paramsList) ) {
    paramsList <- list(paramsList)
  }

  params <- try(modelExtractParam(model, only.values=FALSE), TRUE)
  if (! is.list(params) )
    params <- kernExtractParam(model, only.values=FALSE)
  
  for ( i in seq(along=paramsList) ) {
    if ( is.character(paramsList[[i]]) ) {
      paramInd <- grep(paramsList[[i]], names(params))
      if ( length(paramInd) == 0 )
        warning(paste("No matches for parameter tie spec:", paramsList[[i]]))
    }
    else
      paramInd <- sort(paramsList[[i]])

    if ( any(paramInd[1]==columnToDel) )
      stop("Parameters have already been tied.")

    if ( length(paramInd) > 1 )
      for ( j in seq(2,length.out=(length(paramInd)-1)) ) {
        model$paramGroups[paramInd[j], paramInd[1]] <- 1
        if ( any(paramInd[j]==columnToDel) )
          stop("Parameters have already been tied.")
        columnToDel <- c(columnToDel, paramInd[j])
      }
  }

  if (length(columnToDel) > 0)
    model$paramGroups <- model$paramGroups[,-columnToDel]

  if ( "nParams" %in% names(model) ) {
    model$nParams <- dim(model$paramGroups)[2]
  } else if ( "numParams" %in% names(model) ) {
    model$numParams <- dim(model$paramGroups)[2]
  }

  return (model)
}



.jitCholInv <- function ( M, Num=10, silent=FALSE ) {
  jitter <- 0
  jitter1 <- abs(mean(diag(M)))*1e-6
  eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )

  for ( i in 1:Num ) {

    ## clear the last error message
    try(stop(""),TRUE)

    Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )

    nPos <- grep("not positive definite",  geterrmessage())

    if ( length(nPos) != 0 ) {
      jitter1 <- jitter1*10
      jitter <- jitter1

      if (! silent) {
        warnmsg <- paste("Matrix is not positive definite, adding",
                         signif(jitter,digits=4), "jitter!")
        warning(warnmsg)
      }
    }
    else break
  }

  try(stop(""),TRUE)

  try (solve( Ch, eyeM ), silent=TRUE)
  SingErr <- grep("singular",  geterrmessage())

  if ( length(SingErr) != 0 ) {
    return (NaN)
  }
  else {
    invCh <- solve( Ch, eyeM )
    invM <- invCh %*% t(invCh)

    if ( jitter == 0 ) {
      ans <- list(invM=invM, jitter=jitter, chol=Ch)
    }
    else ans <- list(invM=invM, jitM=M+jitter*eyeM , jitter=jitter, chol=Ch)

    return (ans)
  }
}



.dist2 <- function (x, x2) {
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))

  xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))

  if ( xdim[2] != x2dim[2] )
    stop("Data dimensions are not matched.")

  n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)

  return (n2)
}


#erff <- function(x) 2 * pnorm(x * sqrt(2)) - 1
#erfcf <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
#logerfc <- function(x) log(2) + pnorm(x * sqrt(2), lower = FALSE, log=TRUE)
#erfcx <- function(x) exp(x^2 + logerfc(x))


#miscCDLL <- dyn.load(paste("miscCFunctions", .Platform$dynlib.ext, sep=""))
#ClnDiffErfs <- getNativeSymbolInfo("ClnDiffErfs", miscCDLL)

lnDiffErfs <- function(x1, x2) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  outdim <- pmax(dim(x1), dim(x2))
  outlen <- max(length(x1), length(x2))
  v <- .C("_ClnDiffErfs", as.double(x1), as.double(x2), as.integer(length(x1)), as.integer(length(x2)), v=double(outlen), s=integer(outlen))
  return (list(matrix(data=v$v, outdim), matrix(data=v$s, outdim)));
}


.gradLnDiffErfs <- function(x1, x2, fact1, fact2) {
  m <- pmin(as.matrix(x1)^2, as.matrix(x2)^2)
  dlnPart <- 2/sqrt(pi) * (exp(-x1^2 + m) * fact1 - exp(-x2^2 + m) * fact2)

  g <- list(dlnPart=dlnPart, m=m)
  return (g)

}



modelExtractParam <- function (model, only.values=TRUE) {
  if (is.GPModel(model))
    model <- modelStruct(model)
  
  funcName <- paste(model$type, "ExtractParam", sep="")
  func <- get(funcName, mode="function")
  params <- func(model, only.values)

  if ( !only.values ) {
    origNames <- names(params)
    if ( "paramGroups" %in% names(model) ) {
      paramGroups <- model$paramGroups
      for ( i in seq(length.out=dim(paramGroups)[2]) ) {
        ind <- grep(1, paramGroups[,i])
        if ( is.list(params) ) {
          names(params)[i] <- origNames[ind[1]]
          for ( j in seq(2, length.out=length(ind)-1) )
            names(params)[i] <- paste(names(params)[i], origNames[ind[j]],sep="/")
        }

        paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
      }

      params <- params%*%paramGroups
    }
  }
  return (params)
}



modelExpandParam <- function (model, params) {
  if (is.GPModel(model))
    return (modelExpandParam(modelStruct(model), params))

  if ( is.list(params) )
    params <- params$values

  if ( "paramGroups" %in% names(model) )
    params <- params %*% t(model$paramGroups)

  funcName <- paste(model$type, "ExpandParam", sep="")
  func <- get(funcName, mode="function")
  model <- func(model, params)

  return (model)
}


modelDisplay <- function(model, ...) {
  if (is.GPModel(model))
    model <- modelStruct(model)

  funcName <- paste(model$type, "Display", sep="")
  if(exists(funcName, mode="function")) {
    func <- get(funcName, mode="function")
    func(model, ...)
  }
}

modelObjective <- function (params, model, ...) {
  if (is.GPModel(model))
    model <- modelStruct(model)

  funcName <- paste(model$type, "Objective", sep="")
  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    err <- func(params, model, ...)
  } else {
    funcName <- paste(model$type, "ExpandParam", sep="")
    func <- get(funcName, mode="function")
    model <- func(model, params)

    funcName <- paste(model$type, "LogLikelihood", sep="")
    func <- get(funcName, mode="function")
    err <- - func(model, ...)
  }

  return (err)
}



modelGradient <- function (params, model, ...) {
  if (is.GPModel(model))
    model <- modelStruct(model)

  funcName <- paste(model$type, "Gradient", sep="")

  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    g <- func(params, model, ...)
  } else {
    funcName <- paste(model$type, "ExpandParam", sep="")
    func <- get(funcName, mode="function")
    model <- func(model, params)

    funcName <- paste(model$type, "LogLikeGradients", sep="")
    func <- get(funcName, mode="function")
    g <- - func(model, ...)
  }

  return (g)
}



modelUpdateProcesses <- function (model, predt=NULL) {
  if (is.GPModel(model))
    return (modelUpdateProcesses(modelStruct(model), predt=predt))

  funcName <- paste(model$type, "UpdateProcesses", sep="")
  func <- get(funcName, mode="function")
  return (func(model, predt=predt))
}



modelLogLikelihood <- function (model) {
  if (is.GPModel(model))
    model <- modelStruct(model)

  funcName <- paste(model$type, "LogLikelihood", sep="")
  func <- get(funcName, mode="function")
  return (func(model))
}



##matrixTrace <- function (x) {
##  return ( sum(diag(as.matrix(x))) )
##}



.complexLog <- function (x) {
  if ( is.real(x) & x>0 ) {
      y <- log(x)
  } else {
      if ( is.real(x) & x<0 )
          warning("Log of negative real number, using complex log!")
      y <- log(x+0i)
  }
  return ( y )
}


.distfit <- function(data, dist = "normal") {
  if (dist == "gamma") {
    cdf <- qgamma 
  }

  else if (dist == "normal") {
    cdf <- qnorm
  }

  else {
    stop("Unknown distribution.")
  }

  t <- optim(c(1, 1), fn=.distfit_obj, gr=NULL, data, cdf)

  return (t)
}



.distfit_obj <- function(theta, y, cdf) {
  p <- c(.05, .25, .50, .75, .95)
  x <- cdf(p, theta[1], theta[2])
  r <- .5 * sum((x - y)^2)

  return (r)
}


# Adapted from selectSomeIndex for AnnotatedDataFrame in Biobase
# Original (c) R. Gentleman, V. Carey, M. Morgan, S. Falcon
.listSelectSomeIndex <- function(object, maxToShow=5) {
  len <- length(object)
  if (maxToShow < 3) maxToShow <- 3
  if (len > maxToShow) {
    maxToShow <- maxToShow - 1
    bot <- ceiling(maxToShow/2)
    top <- len-(maxToShow-bot-1)
    list(1:bot, "...", top:len)
  } else if (len >= 1) list(1:len, NULL, NULL)
  else list(NULL, NULL, NULL)
}
