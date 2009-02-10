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
      if ( x[ind] < eps ) y[ind]<- log(eps) else
      y[ind] <- log(x[ind])
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

    if ( any(paramInd==columnToDel) )
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

      if ( any(paramInd==columnToDel) )
        stop("Parameters have already been tied.")

      if ( length(paramsInd) > 0 ) 
        for ( j in seq(2,length.out=(length(paramsInd)-1)) ) {
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
             


erf <- function (x) {
  result <- 2*pnorm(sqrt(2)*x) -1
  return (result)
}



erfc <- function (x) {
  return (1-erf(x))
}



erfcx <- function (x) {
  return ( exp(x*x)*erfc(x) )
}



lnDiffErfs <- function(x1, x2) {
  ## return v = log(erf(x1)-erf(x2))
  x1 <- Re(x1)
  x2 <- Re(x2)
  dim1 <- dim(as.matrix(x1))
  dim2 <- dim(as.matrix(x2))

  v <- matrix(0, max(dim1[1],dim2[1]),  max(dim1[2],dim2[2]))

  if ( length(x1) == 1 )
    x1 <- array(x1, dim=dim2)

  if ( length(x2) == 1 )
    x2 <- array(x2, dim=dim1)

  ## case 1: different signs
  stat1 <- sign(x1)!=sign(x2)
  I1 <- grep(TRUE, stat1)
  ## case 2: both positive signs
  ## 2a: x1>x2
  stat2a <- (x1>0) & (x1>x2) & !stat1
  I2a <- grep(TRUE, stat2a)
  ## 2b: x1<=x2
  stat2b <- (x1>0) & !stat1 & !stat2a
  I2b <- grep(TRUE, stat2b)
  ## case 3: both negative signs
  ## 3a: x1<x2
  stat3a = !stat1 & !stat2a & !stat2b & (x1 < x2);
  I3a <- grep(TRUE, stat3a)
  ## 3b: x1>=x2
  stat3b = !stat1 & !stat2a & !stat2b & !stat3a;
  I3b <- grep(TRUE, stat3b)  

  options(warn=-1)
  if ( length(I1) )
    v[I1] <- expTransform(erf(x1[I1]) - erf(x2[I1]), "xtoa")
  if ( length(I2a) )  
    v[I2a] <- expTransform(erfcx(x2[I2a])-erfcx(x1[I2a])*exp(x2[I2a]^2-x1[I2a]^2), "xtoa")-x2[I2a]^2
  if ( length(I2b) )
    v[I2b] <- expTransform(erfcx(x2[I2b])*exp(x1[I2b]^2-x2[I2b]^2)-erfcx(x1[I2b]), "xtoa")-x1[I2b]^2
  if ( length(I3a) )
    v[I3a] <- expTransform(erfcx(-x1[I3a])*exp(x2[I3a]^2-x1[I3a]^2)-erfcx(-x2[I3a]), "xtoa")-x2[I3a]^2
  if ( length(I3b) )
    v[I3b] <- expTransform(erfcx(-x1[I3b])-erfcx(-x2[I3b])*exp(x1[I3b]^2-x2[I3b]^2), "xtoa")-x1[I3b]^2
  options(warn=0)  

  return(v)

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



modelExpandparam <- function (model, params) {
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
