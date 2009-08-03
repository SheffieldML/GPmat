logLikelihood <- function (model) {
  dataLocation <- model$comp[[1]]
  dim <- length(dataLocation$y)

  ll <- -dim*log(2*pi) - dataLocation$logDetK - t(dataLocation$m) %*% dataLocation$invK %*% dataLocation$m
  ll <- 0.5*ll

  ## prior contributions
  if ( any(grep("bprior",names(model))) ) {
    ll <- ll + kernPriorLogProb(dataLocation$kern)
    ll <- ll + priorLogProb(dataLocation$bprior, dataLocation$B)
  }
  return (ll)
}
