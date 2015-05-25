demAutoOptimiseGp <- function(path=getwd(), filename='demAutoOptimiseGp', png=FALSE, gif=FALSE) {
## DEMAUTOOPTIMISEGP Shows that there is an optimum for the covariance function length-scale.
## DESC Shows that by varying the length scale, an artificial data set has
## different likelihoods, yet there is an optimum for which the likelihood is
## maximised. This demo is similar to demOptimiseGp, only here, it is
## demonstrated how the length scale hyperparameter is optimised automatically
## through SCG (scaled conjugate gradients) numerical optimisation.
## Run multiple times to understand the effect of optimisation on randomly generated datasets.

## COPYRIGHT : Alfredo A. Kalaitzis 2011

## GP
  ## Set up model.
  options = gpOptions()
  options$kern$comp = list("rbf","white")
  lengthScale = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16)

  ## Inputs
  x = matrix(seq(-1, 1, length=6), ncol=1)
  xtest = matrix(seq(-1.5, 1.5, length=200),ncol=1)
  ## Generate y
  trueKern = kernCreate(x, list(type="cmpnd",comp=list("rbf", "white")))
  K = kernCompute(trueKern, x)
  y = t(gaussSamp(matrix(0,6,1), K, 1))
  y = scale(y,scale=FALSE)
  model = gpCreate(dim(x)[2], dim(y)[2], x, y, options)	

  graphics.off();

  ## hyperparameters: inverse-lengthscale, signal-variance, noise-variance.
  inithypers = log(c(1/(lengthScale[1]), 1, model$kern$comp[[2]]$variance))
  model = gpExpandParam(model, inithypers) ## This forces kernel computation.
  ll_init = gpLogLikelihood(model) ## GP log-marginal likelihood for this model.
  meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return=TRUE) ## GP mean and variance.
  dev.new(); plot.new()
  gpPlot(model, xtest, meanVar$mu, meanVar$varsigma, ylim=c(-2.5,2.5), xlim=range(xtest), col='black')
  title(sprintf("l=%.2f, s_f=%.2f, s_n=%.2f\n GP log-likelihood=%.2f",
    1/exp(inithypers[1]), exp(inithypers[2]), exp(inithypers[3]), ll_init));
  if (png) {
    dev.copy2eps(file=paste(path,'/',filename,'_1.eps', sep='') ) ## Save plot as eps
    ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
    system(paste('eps2png ',path,'/',filename,'_1.eps',sep=''))
  }

  model = gpOptimise(model, display=TRUE, iters=400)
  opthypers = gpExtractParam(model, only.values = FALSE);
  opthypers = exp(opthypers);
  ll_opt = gpLogLikelihood(model)
  meanVar = gpPosteriorMeanVar(model, xtest, varsigma.return=TRUE) ## GP mean and variance.
  dev.new(); plot.new()
  gpPlot(model, xtest, meanVar$mu, meanVar$varsigma, ylim=c(-2.5,2.5), xlim=range(xtest), col='black')
  title(sprintf("l=%.2f, s_f=%.2f, s_n=%.2f\n GP log-likelihood=%.2f",
    1/opthypers[1], opthypers[2], opthypers[3], ll_opt));
}
