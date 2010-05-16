gpExtractParam <- function(model, only.values=TRUE) {
  funcName <- optimiDefaultConstraint(model$bTransform)
  func <- get(funcName$func, mode="function")

  if ( only.values ) {
    params <- kernExtractParam(model$kern)
    ## Note: ignores funcName$hasArgs
#     params <- c(params, func(model$B, "xtoa"))
  } else {
    params <- kernExtractParam(model$kern, only.values)
    ## Note: ignores funcName$hasArgs
#     Bparams <- func(model$B, "xtoa")
#     for ( i in seq(along=Bparams) ) {
#       names(Bparams)[i] <- paste("Basal", i, sep="")
#     }
#     params <- c(params, Bparams)
  }

  if ( "fix" %in% names(model) )
    for ( i in seq(along=model$fix$index) )
      params[model$fix$index[i]] <- model$fix$value[i]

  params <- Re(params)

  return (params)
}