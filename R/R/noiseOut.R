## Not implemented: cmpndNoiseOut,mgaussianNoiseOut,ncnmNoiseOut,orderedNoiseOut,
## probitNoiseOut,scaleNoiseOut
noiseOut <- function(noise, mu, varsigma) {
  fhandle = get(paste(noise$type, "NoiseOut", sep=""), mode="function")
  #fhandle = str2func([noise.type 'NoiseOut']);
  y = fhandle(noise, mu, varsigma)

  return (y)
}
