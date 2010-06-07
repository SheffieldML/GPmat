gpPlot <- function(model, Xstar, mu, S, simpose=NULL, col='blue', title) {
  
  if (missing(model) || missing(Xstar)) {
    stop('Missing GP model or points of prediction Xstar.')
  } else {
    if (missing(mu) || missing(S)) {
      meanVar = gpPosteriorMeanVar(model, Xstar, varsigma.return=TRUE)
      mu = meanVar$mu; S = meanVar$varsigma
    }
  }

#   par(pty="s") ## set plot basis
  plot(0, type="n", xlim=c(-10,250), ylim=c(-0.5, 0.5), xlab='', ylab='',main=title)

  f = c(mu+2*sqrt(S), rev(mu-2*sqrt(S)))
  polygon(c(Xstar, rev(Xstar)), f, col = "grey", border = NA)
  lines(Xstar, mu, col=col, lwd=.3)
  points(model$X, model$y, pch = 3, cex = 1.5, col = col)

  if (!is.null(simpose)) {
    y = mu[simpose] + rnorm(6, 0, exp(model$params$xmin[3]/2))
    points(simpose, y, pch = 4, cex = 1.5, lwd=3, col = col)
  }
}