kernDiagGradient <- function (kern, x, covDiag) {
  x = as.array(x)
  funcName <- paste(kern$type, 'KernDiagGradient', sep="")

  if (exists(funcName)) {
    fhandle = get(funcName, mode="function")
    g = fhandle(kern, x, covDiag)
  } else {
    fhandle = get(paste(kern$type, 'KernGradient', sep=""))
    g = matrix(0, 1, kern$nParams)
    for (i in 1:tail(dim(x),1)) {
      g = g + fhandle(kern, x[i, ], covDiag[i])
    }
  }

  ## Check if parameters are being optimised in a transformed space.
  factors = .kernFactors(kern, 'gradfact')

  for (i in seq(along=factors))
    g[factors[[i]]$index] <- g[factors[[i]]$index]*factors[[i]]$val

  return (g)
}