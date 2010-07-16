priorCreate <- function(type) {
  if (!is.list(type)) {
    prior <- list(type=type)
    return (priorParamInit(prior))
  } else {
    prior <- list(type=type$type)
    prior <- priorParamInit(prior)
    if ("index" %in% names(type))
      prior$index <- type$index
    return (priorExpandParam(prior, type$params, untransformed.values=TRUE))
  }
}


priorParamInit <- function(prior) {
  func <- get(paste(prior$type, 'PriorParamInit', sep=''), mode='function')
  return (func(prior))
}


priorExpandParam <- function (prior, params, untransformed.values=FALSE) {
  if ( is.list(params) )
    params <- params$values
  
  if ( "transforms" %in% names(prior) && (length(prior$transforms) > 0)
      && !untransformed.values )
    for ( i in seq(along=prior$transforms) ) {
      index <- prior$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(prior$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "atox", prior$transformArgs[[i]])
      else
        params[index] <- func(params[index], "atox")
    }

  funcName <- paste(prior$type, "PriorExpandParam", sep="")
  func <- get(funcName, mode="function")
  return (func(prior, params))
}


priorExtractParam <- function (prior, only.values=TRUE, untransformed.values=FALSE) {
  funcName <- paste(prior$type, "PriorExtractParam", sep="")
  func <- get(funcName, mode="function")

  params <- func(prior, only.values=only.values, untransformed.values=untransformed.values)

  if ( any(is.nan(params)) )
    warning("Parameter has gone to NaN.")

  if ( "transforms" %in% names(prior) && (length(prior$transforms) > 0)
      && !untransformed.values )
    for ( i in seq(along=prior$transforms) ) {
      index <- prior$transforms[[i]]$index
      funcName <- optimiDefaultConstraint(prior$transforms[[i]]$type)
      func <- get(funcName$func, mode="function")
      if (funcName$hasArgs)
        params[index] <- func(params[index], "xtoa", prior$transformArgs[[i]])
      else
        params[index] <- func(params[index], "xtoa")
    }

  return (params)
}


priorGradient <- function(prior, params) {
  func <- get(paste(prior$type, 'PriorGradient', sep=''), mode='function')
  return (func(prior, params))
}


priorLogProb <- function(prior, x) {
  func <- get(paste(prior$type, 'PriorLogProb', sep=''), mode='function')
  return (func(prior, x))
}
