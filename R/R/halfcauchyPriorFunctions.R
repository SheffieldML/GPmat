halfcauchyPriorExpandParam <- function(prior, params) {
  prior$s <- params[1]

  return (prior)
}


halfcauchyPriorExtractParam <- function(prior, only.values=TRUE,
                                   untransformed.values=TRUE) {
  params <- c(prior$s)

  if ( !only.values )
    names(params) <- c("half-cauchy scale")

  return (params)
}


halfcauchyPriorGradient <- function(prior, x) {
  return (- 4 * prior$s * x / (pi * (x*x + (prior$s)^2)))
}


halfcauchyPriorLogProb <- function(prior, x) {
  D <- length(x)
  return (D*log(2*prior$s) - D*log(pi) - sum(log(x*x + prior$s^2)))
}


halfcauchyPriorParamInit <- function(prior) {
  prior$s <- 1

  prior$transforms <- list(list(index=c(1), type="positive"))
  prior$nParams <- 1

  return (prior)
}
