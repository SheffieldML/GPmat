gpPlot <- function(model,Xstar,mu,S,simpose=NULL,xlim=NULL,ylim=NULL,xlab='',ylab='',col='blue',title='') {
## GPPLOT Plots the GP mean and variance.

## COPYRIGHT : Alfredo A. Kalaitzis 2010

## GP

  if (missing(model) || missing(Xstar)) {
    stop('Missing GP model or points of prediction Xstar.')
  } else {
    if (missing(mu) || missing(S)) {
      meanVar = gpPosteriorMeanVar(model, Xstar, varsigma.return=TRUE)
      mu = meanVar$mu; S = meanVar$varsigma
    }
  }

#   f = c(mu+2*sqrt(abs(S)), rev(mu-2*sqrt(abs(S))))
  f = c(mu+2*abs(S), rev(mu-2*abs(S)))

  xcomb = rbind(model$X,Xstar); ycomb = f # rbind(model$y,mu,S)
  if (is.null(xlim))
    xlim = range(xcomb) #c(min(xcomb), max(xcomb))
  if (is.null(ylim))
    ylim = range(ycomb) #c(min(ycomb), max(ycomb))

#   par(pty="s") 
  plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=title, new=TRUE) ## Empty plot basis.

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
