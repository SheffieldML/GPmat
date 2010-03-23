

.baselineOptimise <- function(y, yvar, options) {
  if (options$includeNoise) {
    oldsigma <- 10
    sigma <- 0
    while (abs(sigma - oldsigma) > 1e-10) {
      oldsigma <- sigma
      grad <- .baselineLogLikeGradient(expTransform(sigma, 'atox'), y, yvar)
      sigma <- sigma - expTransform(sigma, 'atox') * grad$d1 /
        (expTransform(2*sigma, 'atox') * grad$d2
         + expTransform(sigma, 'atox') * grad$d1)
    }
    sigma <- expTransform(sigma, 'atox')
  }
  else
    sigma <- 0

  return (.baselineLogLikelihood(sigma, y, yvar))
}


.baselineLogLikeGradient <- function(sigma, y, yvar) {
  mu <- sum(y / (yvar + sigma)) / sum(1/(yvar+sigma))
  return (list(d1=0.5 * sum((y - mu)^2 / (yvar + sigma)^2 - 1/(yvar + sigma)),
               d2=0.5 * sum(1/(yvar + sigma)^2) - sum((y - mu)^2 / (yvar + sigma)^3)))
}


.baselineLogLikelihood <- function(sigma, y, yvar) {
  mu <- sum(y / (yvar + sigma)) / sum(1/(yvar+sigma))
  return (-0.5 * sum(log(2*pi*(yvar+sigma)) + (y - mu)^2 / (yvar + sigma)))
}
