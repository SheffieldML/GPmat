distfit <- function(data, dist = "normal") {

  if (dist == "gamma") {
    cdf <- qgamma 
  }

  else if (dist == "normal") {
    cdf <- qnorm
  }

  else {
    stop("Unknown distribution.")
  }

  t <- optim(c(1, 1), fn=distfit_obj, gr=NULL, data, cdf)

  return (t)
}


distfit_obj <- function(theta, y, cdf) {

  p <- c(.05, .25, .50, .75, .95)
  x <- cdf(p, theta[1], theta[2])
  r <- .5 * sum((x - y)^2)

  return (r)
}
