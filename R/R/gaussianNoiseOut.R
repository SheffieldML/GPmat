gaussianNoiseOut <- function(noise, mu, varsigma) {
  D = dim(mu)[2]
  y = mu * 0
  for (i in 1:D)
    y[,i] = mu[,i] + noise$bias[i]

  return (y)
}
