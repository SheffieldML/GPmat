

.baselineOptimise <- function(y, yvar, options) {
  if (options$includeNoise) {
    # Check for special case of constant expression
    if (sd(y) == 0) {
      return (Inf);
    }
    oldsigma <- 10
    newsigma <- 0
    oldval <- .baselineLogLikelihood(expTransform(newsigma, 'atox'), y, yvar)
    l <- 0
    while ((abs(newsigma - oldsigma) > 1e-10) && l < 1000) {
      oldsigma <- newsigma
      grad <- .baselineLogLikeGradient(expTransform(oldsigma, 'atox'), y, yvar)
      step <- expTransform(oldsigma, 'atox') * grad$d1 /
        (expTransform(2*oldsigma, 'atox') * grad$d2
         + expTransform(oldsigma, 'atox') * grad$d1)
      k <- 0
      newval <- oldval - 10
      while (newval < (oldval - 1e-10) && k < 20) {
        newsigma <- oldsigma - step
        newval <- .baselineLogLikelihood(expTransform(newsigma, 'atox'), y, yvar)
        step <- step / 2
        k <- k+1
      }
      oldval <- newval
      l <- l+1
    }
    sigma <- expTransform(newsigma, 'atox')
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
