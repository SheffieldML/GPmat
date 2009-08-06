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

  if ( is.list(kernType) & any(grep("complete",names(kernType))) ) {
    if ( kernType$complete == 1 ) {
      kern <- kernType
    }
    
  } else if ( is.list(kernType) ) {
    
    kern <- list(inputDimension=dim, type=kernType$type)
    
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
        
    } else if ( any(grep(kern$type, c("cmpnd", "tensor", "translate"))) )  {
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
      if ( start == length(kernType) ) {
        kern$argument <- kernCreate(x, kernType$comp[start])
      } else {
        kern$argument <- kernCreate(x, kernType$comp[start:length(kernType)])
      }
    }

    kern <- kernParamInit(kern)

  } else {
    kern <- list(type=kernType, inputDimension=dim)
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

    if ( any(grep("transforms", names(kern))) ) 
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

    if ( any(grep("transforms", names(kern))) ) 
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
  
  if ( any(grep("transforms", names(kern))) ) 
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
