invgammaPriorExpandParam <- function(prior, params) {
  prior$a <- params[1]
  prior$b <- params[2]

  return (prior)
}


invgammaPriorExtractParam <- function(prior, only.values=TRUE,
                                   untransformed.values=TRUE) {
  params <- c(prior$a, prior$b)

  if ( !only.values )
    names(params) <- c("inv gamma a", "inv gamma b")

  return (params)
}


invgammaPriorGradient <- function(prior, x) {
  return (- (prior$a + 1) / x + prior$b / (x*x))
}


invgammaPriorLogProb <- function(prior, x) {
  D <- length(x)
  return (D*prior$a*log(prior$b) - D*lgamma(prior$a) - (prior$a+1)*sum(log(x))
          - prior$b*sum(1/x))
}


invgammaPriorParamInit <- function(prior) {
  prior$a <- 1e-6
  prior$b <- 1e-6

  prior$transforms <- list(list(index=c(1,2), type="positive"))
  prior$nParams <- 2

  return (prior)
}
