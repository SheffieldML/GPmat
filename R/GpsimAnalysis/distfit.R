distfit <- function(data, dist = "normal") {

  source("distfit_obj.R")

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
