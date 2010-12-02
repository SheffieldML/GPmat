noiseCreate <- function(noiseType, y) {
  if (is.array(noiseType)) #isstruct(noiseType)
    return (noiseType)
  else if (is.list(noiseType)) { #iscell(noiseType)
    ## compound noise type
    noise$type = "cmpnd"
    if (nargs() > 1) {
      for (i in 1:length(noiseType))
	noise$comp[[i]] = noiseCreate(noiseType[[i]], y[,i])
    } else {
      for (i in 1:length(noiseType))
	noise$comp[[i]] = noiseCreate(noiseType[[i]])
    }
  } else {
    noise$type = noiseType
  }

  if (nargs() > 1)
    noise = noiseParamInit(noise, y)

  ## Check if the noise model has bespoke site update code
  noise$updateSites = 0
  if (exists(paste(noise$type,"NoiseSites",sep="")))
    noise$updateSites = 1

  ## Check if the model has bespoke nu and g update code.
  noise$updateNuG = 0
  if (exist(paste(noise$type,"NoiseNuG",sep="")))
    noise$updateNuG = 1

  return (noise)
}