expTransform <- function (x, transform="atox") {

  eps <- 2.2204e-16
  maxVal <- 4.3112e+15

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( any(grep("atox", transform)) ) {
    for ( ind in 1:length(as.array(x)) ) {
      if ( x[ind] > thre ) y[ind] <- maxVal else
      if ( x[ind] < -thre ) y[ind]<- eps else
      y[ind] <- exp(x[ind])
    }
  } else
  if ( any(grep("xtoa", transform)) ) {
    for ( ind in 1:length(as.array(x)) ) {
        y[ind] <- complexLog(x[ind])
    }
  } else
  if ( any(grep("gradfact", transform)) )
    y <- x

  return (y)
}


modelTieParam <- function (model, paramsList) {
  columnToDel <- c()

  if ( !is.list(paramsList) & length(paramsList)>0 ) {

    paramInd <- sort(paramsList)

    if ( any(paramInd[1]==columnToDel) )
      stop("Parameters have already been tied.")

    for ( j in seq(2,length.out=(length(paramInd)-1)) ) {
      model$paramGroups[paramInd[j], paramInd[1]] <- 1
      if ( any(paramInd[j]==columnToDel) )
        stop("Parameters have already been tied.")
      columnToDel <- c(columnToDel, paramInd[j])
    }

  } else if ( is.list(paramsList) ) {

    for ( i in seq(length.out=length(paramsList)) ) {

      paramInd <- sort(paramsList[[i]])

      if ( any(paramInd[1]==columnToDel) )
        stop("Parameters have already been tied.")

      if ( length(paramInd) > 0 )
        for ( j in seq(2,length.out=(length(paramInd)-1)) ) {
          model$paramGroups[paramInd[j], paramInd[1]] <- 1
          if ( any(paramInd[j]==columnToDel) )
            stop("Parameters have already been tied.")
          columnToDel <- c(columnToDel, paramInd[j])
        }
    }
  }

  model$paramGroups <- model$paramGroups[,-columnToDel]

  if ( any(grep("nParams", names(model))) ) {
    model$nParams <- dim(model$paramGroups)[2]
  } else if ( any(grep("numParams", names(model))) ) {
    model$numParams <- dim(model$paramGroups)[2]
  }

  return (model)
}



jitCholInv <- function ( M, Num=10 ) {
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

      warnmsg <- paste("Matrix is not positive definite, adding",
                       signif(jitter,digits=4), "jitter!")
      warning(warnmsg)
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



dist2 <- function (x, x2) {
  xdim <- dim(as.matrix(x))
  x2dim <- dim(as.matrix(x2))

  xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
  x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))

  if ( xdim[2] != x2dim[2] )
    stop("Data dimensions are not matched.")

  n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)

  return (n2)
}


erff <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfcf <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
logerfc <- function(x) log(2) + pnorm(x * sqrt(2), lower = FALSE, log=TRUE)
erfcx <- function(x) exp(x^2 + logerfc(x))


#miscCDLL <- dyn.load(paste("miscCFunctions", .Platform$dynlib.ext, sep=""))
#ClnDiffErfs <- getNativeSymbolInfo("ClnDiffErfs", miscCDLL)

lnDiffErfs <- function(x1, x2) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  outdim <- pmax(dim(x1), dim(x2))
  outlen <- max(length(x1), length(x2))
  v <- .C(ClnDiffErfs, as.double(x1), as.double(x2), as.integer(length(x1)), as.integer(length(x2)), v=double(outlen), s=integer(outlen))
  return (list(matrix(data=v$v, outdim), matrix(data=v$s, outdim)));
}


lnDiffErfs2 <- function(x1, x2) {
 
  ## return v = log(erff(x1)-erff(x2))
  x1 <- Re(x1)
  x2 <- Re(x2)
  dim1 <- dim(as.matrix(x1))
  dim2 <- dim(as.matrix(x2))

  v <- matrix(0, max(dim1[1],dim2[1]),  max(dim1[2],dim2[2]))

  if ( length(x1) == 1 )
    x1 <- array(x1, dim=dim2)

  if ( length(x2) == 1 )
    x2 <- array(x2, dim=dim1)

  signs = sign(x1 - x2)
  I = (signs == -1)
  swap <- x1[I]
  x1[I] <- x2[I]
  x2[I] <- swap

  ## case 1: different signs
  I1 <- (x1 * x2) < 0
  ## case 2: equal
  I2 <- x1 == x2
  ## case 3: both positive
  I3 <- (x1>0) & !I1 & !I2
  ## case 4: both negative
  I4 = !I1 & !I2 & !I3

  options(warn=-1)
  # The absolute values are to maintain numerical stability, they should only
  # ever be applied for really small arguments if there are numerical problems
  v[I1] <- log(erff(x1[I1]) - erff(x2[I1]))
  v[I2] <- -Inf
  v[I3] <- log(abs(erfcx(x2[I3])-erfcx(x1[I3])*exp(x2[I3]^2-x1[I3]^2)))-x2[I3]^2
  v[I4] <- log(abs(erfcx(-x1[I4])-erfcx(-x2[I4])*exp(x1[I4]^2-x2[I4]^2)))-x1[I4]^2
  options(warn=0)

  return(list(v, signs))
}



gradLnDiffErfs <- function(x1, x2, fact1, fact2) {
  m <- pmin(as.matrix(x1)^2, as.matrix(x2)^2)
  dlnPart <- 2/sqrt(pi) * (exp(-x1^2 + m) * fact1 - exp(-x2^2 + m) * fact2)

  g <- list(dlnPart=dlnPart, m=m)
  return (g)

}



listStruct <- function (list) {
  if ( is.list(list) ) {
    length <- length(list)
    names <- names(list)

    for ( i in seq(length.out=length) ) {
      if ( is.character(list[[i]]) | is.logical(list[[i]]) ) {
        if (length(list[[i]])<10) {
          cat (names[i], ":  ", list[[i]],"\n")
        } else {
          cat (names[i], ":  1x", length(list[[i]]), mode(list[[i]]), "\n")
        }
      } else if ( is.numeric(list[[i]]) ) {
        if ( is.matrix(list[[i]]) ) {
          dim <- dim(list[[i]])
          msg <- paste(names[i], ":  ", dim[1], "x", dim[2], "  ", mode(list[[i]]), sep="")
          cat (msg,"\n")
        } else
          if ( length(list[[i]])<7 ) {
            cat (names[i], ":  [", list[[i]],"] \n")
          } else {
            cat (names[i],":  1x",length(list[[i]]),"  ",mode(list[[i]]), "\n")
          }
      } else if ( is.list(list[[i]]) ) {
        msg <- paste(names[i],":  ","1x",length(list[[i]]),"  ",mode(list[[i]]),sep="")
        cat (msg,"\n")
      } else {
        cat (names[i], ":", mode(list[[i]]), " \n")
      }
    }
  } else warning("the argument is not a list.")
}



modelExtractParam <- function (model, option=1) {
  funcName <- paste(model$type, "ExtractParam", sep="")
  func <- get(funcName, mode="function")
  params <- func(model)

  if ( option>1 ) {
    if ( any(grep("paramGroups", names(model))) ) {
      paramGroups <- model$paramGroups
      for ( i in seq(length.out=dim(paramGroups)[2]) ) {
        ind <- grep(1, paramGroups[,i])
        if ( is.list(params) ) {
          params$names[i] <- origNames[ind[1]]
          for ( j in seq(2, length.out=length(ind)-1) )
            params$names[i] <- paste(params$names[i], origNames[ind[j]],sep="/")
        }

        paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
      }

      if ( is.list(params) ) {
        params$values <- params$values%*%paramGroups
      } else {
        params <- params%*%paramGroups
      }
    }
  }
  return (params)
}



modelExpandParam <- function (model, params) {
  if ( is.list(params) )
    params <- params$values

  if ( any(grep("paramGroups", names(model))) )
    params <- params %*% t(model$paramGroups)

  funcName <- paste(model$type, "ExpandParam", sep="")
  func <- get(funcName, mode="function")
  model <- func(model, params)

  return (model)
}



modelObjective <- function (params, model, ...) {
  funcName <- paste(model$type, "Objective", sep="")
  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    err <- func(params, model, ...)
  } else {
    funcName <- paste(model$type, "ExpandParam", sep="")
    func <- get(funcName, mode="function")
    model <- func(model, params)

    funcName <- paste(model$type, "Likelihood", sep="")
    func <- get(funcName, mode="function")
    err <- - func(model, ...)
  }

  return (err)
}



modelGradient <- function (params, model, ...) {
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



sigmoid <- function (x) {
  y <- array(1, dim(x))/(1+exp(-x))
  return (y)
}



trace <- function (x) {
  return ( sum(diag(as.matrix(x))) )
}



complexLog <- function (x) {
  if ( is.real(x) & x>0 ) {
      y <- log(x)
  } else {
      if ( is.real(x) & x<0 )
          warning("Log of negative real number, using complex log!")
      y <- log(x+0i)
  }
  return ( y )
}


logLikelihood <- function (model) {
  dataLocation <- model$comp[[1]]
  dim <- length(dataLocation$y)

  ll <- -dim*log(2*pi) - dataLocation$logDetK - t(dataLocation$m) %*% dataLocation$invK %*% dataLocation$m
  ll <- 0.5*ll

  ## prior contributions
  if ( any(grep("bprior",names(model))) ) {
    ll <- ll + kernPriorLogProb(dataLocation$kern)
    ll <- ll + priorLogProb(dataLocation$bprior, dataLocation$B)
  }
  return (ll)
}



distfit <- function(data, dist = "normal") {

  if (dist == "gamma") {
    cdf <- qgamma 
  }

  else if (dist == "normal") {
    cdf <- qnorm
  }

  else {
    stop("Unknown distribution.")
  }

  t <- optim(c(1, 1), fn=distfit_obj, gr=NULL, data, cdf)

  return (t)
}



distfit_obj <- function(theta, y, cdf) {

  p <- c(.05, .25, .50, .75, .95)
  x <- cdf(p, theta[1], theta[2])
  r <- .5 * sum((x - y)^2)

  return (r)
}
