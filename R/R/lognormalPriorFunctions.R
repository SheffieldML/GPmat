lognormalPriorExpandParam <- function(prior, params) {
  prior$mu <- params[1]
  prior$s2 <- params[2]

  return (prior)
}


lognormalPriorExtractParam <- function(prior, only.values=TRUE,
                                       untransformed.values=TRUE) {
  params <- c(prior$mu, prior$s2)

  if ( !only.values )
    names(params) <- c("log-normal mean", "log-normal variance")

  return (params)
}


lognormalPriorGradient <- function(prior, x) {
  return (- (log(x)-prior$mu)/(prior$s2 * x))
}


lognormalPriorLogProb <- function(prior, x) {
  D <- length(x)
  return (-0.5*D*log(2*pi*prior$s2) - 0.5*sum((log(x)-prior$mu)^2)/prior$s2)
}


lognormalPriorParamInit <- function(prior) {
  prior$mu <- 0
  prior$s2 <- 1

  prior$transforms <- list(list(index=2, type="positive"))
  prior$nParams <- 2

  return (prior)
}
