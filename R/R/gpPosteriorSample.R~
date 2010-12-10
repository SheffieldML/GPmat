gpPosteriorSample <- function(kernType, numSamps=10, params=NULL, lims=c(-3,3), bw=FALSE, path=getwd()) {
## GPPOSTERIORSAMPLE Create a plot of samples from a posterior covariance.
## FORMAT
## DESC creates a plot of samples from a kernel with the given parameters and variance.
## ARG kernType : the type of kernel to sample from.
## ARG numSamps : the number of samples to take.
## ARG params : parameter vector for the kernel.
## ARG lims : limits of the x axis.
##
## COPYRIGHT : Neil D. Lawrence, 2008; Alfredo A. Kalaitzis, 2010

  t_star = matrix(seq(lims[1], lims[2], length=200), ncol=1)

  kern = kernCreate(t_star, kernType)
  for (i in 1:length(kern$comp)) {
    kern$comp[[i]]$transforms = NULL
  }

  if (!is.null(params)) {
    feval = get(paste(kern$type,'KernExpandParam',sep=''), mode="function")
    kern = feval(kern, params)
  }
  feval = get(paste(kernType,'KernExtractParam',sep=''), mode="function")
  params = feval(kern, only.values=FALSE)
  paramStr = c()

  for (i in 1:length(params)) {
    Name = strsplit(names(params)[i], split='')[[1]]
    Name[1] = toupper(Name[1])
    Name = paste(Name, collapse='')
    ind = grep(' ', strsplit(Name, split='')[[1]])
    Name[ind+1] = toupper(Name[ind+1])
    Name[ind] = ''
    paramStr = paste(paramStr, Name, as.character(params[i]), sep='')
  }
#   paramStr(find(paramStr==46)) = 'p';
  infoStr = paste('Samples', as.character(numSamps), sep='')

  ## Covariance of the prior.
  K_starStar = kernCompute(kern, t_star, t_star)

  ## Generate "training data" from a sine wave.
  x = matrix(runif(5)*(lims[2]-lims[1])*0.6 + (lims[1]*0.6), ncol=1)
  f = sin(x)

  ## Compute kernel for training data.
  K_starf = kernCompute(kern, t_star, x)
  K_ff = kernCompute(kern, x)

  ## Mean and covariance of posterior.
  fbar = K_starf %*% .jitCholInv(K_ff, silent=TRUE)$invM %*% f
  Sigma = K_starStar - (K_starf %*% .jitCholInv(K_ff, silent=TRUE)$invM %*% t(K_starf))

  ## Sample from the posterior.
  fsamp = gaussSamp(fbar, Sigma, numSamps)

  ## Plot and save
  graphics.off() ## kill all devices
  dev.new(width=5,height=4); plot.new()

  matplot(t_star, t(fsamp), ylim=range(fsamp),type='l', xlab='x',ylab='f')
  points(x, f)
  zeroAxes()

  if (is.list(kernType)) {
    KernType = NULL
    for (i in seq(length(kernType),1,-1)) {
      KernType = paste(kernType[[i]], KernType)
      KernType = strsplit(KernType, split='')[[1]]
      KernType[1] = toupper(KernType[1])
      KernType = paste(KernType, collapse='')
    }
  }
  else {
    KernType = strsplit(kernType, split='')[[1]]
    KernType[1] = toupper(KernType[1])
    KernType = paste(KernType, collapse='')
  }

  fileName = paste('gpPosteriorSample', KernType, infoStr, paramStr, sep='')
  dev.copy2eps(file = paste(path,fileName, '.eps', sep='')) ## Save plot as eps.
  ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
  system(paste('eps2png ',path,fileName, '.eps',sep=''))
# 

}