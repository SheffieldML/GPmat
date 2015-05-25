cmpndNoiseParamInit <- function(noise, y) {
  if (nargs() > 1) {
    if (length(noise$comp) != dim(y)[2])
      stop("Number of noise components must match dimensions of y.")
  }
  noise$nParams = 0
  for (i in 1:length(noise$comp)) {
    if (nargs() > 1)
      noise$comp[[i]] = noiseParamInit(noise$comp[[i]], y[,i])
    else {
      noise$comp[[i]]$numProcess = 1
      noise$comp[[i]] = noiseParamInit(noise$comp[[i]])
    }
    noise$nParams = noise$nParams + noise$comp[[i]]$nParams
  }
  noise$paramGroups = diag(1, nrow=noise$nParams, ncol=noise$nParams) #speye(noise.nParams)
  ## This is a bit of a hack.
  noise$missing=0

  return (noise)
}
