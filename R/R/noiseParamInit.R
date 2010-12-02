## Not implemented: mgaussianNoiseParamInit, ncnmNoiseParamInit,
## ngaussNoiseParamInit, noiseParamInit, orderedNoiseParamInit,
## probitNoiseParamInit, scaleNoiseParamInit
noiseParamInit <- function(noise, y) {
  noise$spherical = 0
  noise$logconcave = 1

  fhandle <- get(paste(noise$type, "NoiseParamInit", sep=""), mode="function")
  if (nargs() > 1)
    noise = fhandle(noise, y)
  else
    noise = fhandle(noise)

  return (noise)
}
