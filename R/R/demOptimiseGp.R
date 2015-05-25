demOptimiseGp <- function(path=getwd(), filename='demOptimiseGp', png=FALSE, gif=FALSE) {
## DEMOPTIMISEGP Shows that there is an optimum for the covariance function length scale.
## DESC shows that by varying the length scale an artificial data
## set has different likelihoods, yet there is an optimum for which
## the likelihood is maximised.

## COPYRIGHT : Neil D. Lawrence 2006, 2008, Alfredo A. Kalaitzis 2010

## GP

  ## Set up model.
  options = gpOptions()
  options$kern$comp = list("rbf","white")
  lengthScale = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)
  figNo = 0

  ## Inputs
  x = matrix(seq(-1, 1, length=6), ncol=1)
  xtest = matrix(seq(-1.5, 1.5, length=200),ncol=1)
  ## Generate y
  trueKern = kernCreate(x, list(type="cmpnd",comp=list("rbf", "white")))
  kern = trueKern
  K = kernCompute(trueKern, x)
  y = t(gaussSamp(matrix(0,6,1), K, 1))
  y = scale(y,scale=FALSE)
  model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)	

  graphics.off();

  ll=c(); llLogDet=c(); llFit=c();
  for (i in 1:length(lengthScale)) {
    # kern$comp[[1]]$inverseWidth = 1/(lengthScale[i]^2)
    # K = kernCompute(kern, x)
    # invK = .jitCholInv(K, silent=TRUE)
    # logDetK = 2* sum( log ( diag(invK$chol) ) )
    ## Hyperparameters: inverse-lengthscale, signal-variance, noise-variance.
    inithypers = log(c(1/(lengthScale[i]), 1, model$kern$comp[[2]]$variance))
    model = gpExpandParam(model, inithypers) ## This forces kernel computation. 
    invK = model$invK_uu
    logDetK = model$logDetK_uu
    ##
    # ll[i] = -0.5*(logDetK + t(y)%*%invK%*%y + dim(y)[1]*log(2*pi))
    ll[i] = gpLogLikelihood(model) ## GP log-marginal likelihood for this model.
    llLogDet[i] = -.5 * (logDetK + dim(y)[1] * log(2*pi))
    llFit[i] = -.5 * t(y) %*% invK %*% y
    ##
    # Kx = kernCompute(kern, x, xtest)
    # mu = t(Kx)%*%invK%*%y
    # S = kernDiagCompute(kern, xtest) - rowSums((t(Kx)%*%invK) * t(Kx))
    meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return=TRUE) ## GP mean and variance.
    # mu = meanVar$mu
    # S = meanVar$varsigma
    # S = S - exp(2*inithypers[3]) ## subtract noise variance

    dev.new(); plot.new() #dev.set(2)
    gpPlot(model, xtest, meanVar$mu, meanVar$varsigma, ylim=c(-2.5,2.5), xlim=range(xtest), col='black')
    figNo = figNo + 1
    if (png) {
      dev.copy2eps(file=paste(path,'/',filename,'1_', as.character(figNo), '.eps', sep='') ) ## Save plot as eps
      ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
      system(paste('eps2png ',path,'/',filename,'1_',as.character(figNo),'.eps',sep=''))
    }

    dev.new(); plot.new() #dev.set(3)
    matplot(lengthScale[1:i], cbind(ll[1:i], llLogDet[1:i],  llFit[1:i]), type="l", lty=c(1,2,3), log="x",
      xlab='length-scale', ylab='log-probability')
    legend(x='topleft',c('marginal likelihood','minus complexity penalty','data-fit'),
      lty=c(1,2,3),col=c('black','red','green'))
    if (png) {
      dev.copy2eps(file = paste(path,'/',filename,'2_', as.character(figNo), '.eps', sep='')) ## Save plot as eps
      system(paste('eps2png ',path,'/',filename,'2_', as.character(figNo), '.eps',sep=''))
    }
  }

  ## Convert the .png files to one .gif file using ImageMagick. The -delay flag sets the time between showing
  ## the frames, i.e. the speed of the animation.
  if (gif) {
    system(paste('convert -delay 80 ',path,'/',filename,'1_*.png ', path,'/',filename,'1.gif', sep=''))
    system(paste('convert -delay 80 ',path,'/',filename,'2_*.png ', path,'/',filename,'2.gif', sep=''))
  }
}
