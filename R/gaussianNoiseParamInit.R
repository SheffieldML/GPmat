gaussianNoiseParamInit <- function(noise, y) {
  if (nargs() > 1) {
    noise$bias = mean(y)
    noise$numProcess = dim(y)[2]
  } else {
    noise$bias = matrix(0, 1, noise$numProcess)
  }

  noise$sigma2 = 1e-6

  noise$transforms <- list(list(index=c(noise$numProcess+1), type="positive"))
  #noise.transforms.index = noise.numProcess+1;
  noise.nParams = 1 + noise$numProcess

  ## Can handle missing values?
  noise$missing = 0;

  ## Noise model leads to constant value of beta.
  noise$spherical = 1

  return (noise)
}