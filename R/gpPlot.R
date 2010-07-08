gpPlot <- function(model,Xstar,mu,S,simpose=NULL,xlim=NULL,ylim=NULL,col='blue',title='') {

  if (missing(model) || missing(Xstar)) {
    stop('Missing GP model or points of prediction Xstar.')
  } else {
    if (missing(mu) || missing(S)) {
      meanVar = gpPosteriorMeanVar(model, Xstar, varsigma.return=TRUE)
      mu = meanVar$mu; S = meanVar$varsigma
    }
  }

  if (is.null(xlim))
    xlim=c(min(rbind(model$X,Xstar)), max(rbind(model$X,Xstar)))
  if (is.null(ylim))
    ylim=c(min(rbind(model$y,mu)), max(rbind(model$y,mu)))

#   par(pty="s") ## set plot basis
  plot(0, type="n", xlim=xlim, ylim=ylim, xlab='', ylab='',main=title)

  f = c(mu+2*sqrt(abs(S)), rev(mu-2*sqrt(abs(S))))
  if (col=='blue') shade = rgb(0,0,255,max=255,alpha=.1*255)
  else if (col=='red') shade = rgb(255,0,0,max=255,alpha=.1*255)
  else shade = 'gray'
  polygon(c(Xstar, rev(Xstar)), f, col = shade, border = shade)
  lines(Xstar, mu, col=col, lwd=2)
  points(model$X, model$y, pch = 3, cex = .5, lwd=2, col = col)

  if (!is.null(simpose)) {
    y = mu[simpose] + rnorm(6, 0, exp(model$params$xmin[3]/2))
    points(simpose, y, pch = 4, cex = 1.5, lwd=3, col = col)
  }

  zeroAxes(cbind(model$X,model$y))
}
