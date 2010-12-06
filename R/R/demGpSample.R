## DEMGPSAMPLE Simple demonstration of sampling from a covariance function.
demGpSample <- function( bw=FALSE, path=getwd(), filename='gpSample') {
  require(fields)

  x = as.matrix(seq(-1,1,length=25))
  kern = kernCreate(x, 'rbf')
  kern$inverseWidth = 10

  graphics.off() ## kill all devices
  par(pty="s")
  par(mar=c(4, 4, 2, 6)) ## increase right margin for legend
  K = kernCompute(kern, x)
  colormap=heat.colors(256)
  if (bw) colormap=gray.colors(256)
  
  image(1:dim(K)[2], 1:dim(K)[1] , K[dim(K)[2]:1, ], xlab = "n", ylab = "n")
  image.plot(legend.only = TRUE, nlevel = 256, zlim = c(min(K),max(K)), col = colormap) ## Colorbar
  dev.copy2eps(file = paste(path,'gpCovariance.eps',sep='')) ## Save plot as eps.
  ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
  system(paste('eps2png ', path, 'gpCovariance.eps', sep=''))

  ## Need to take the real part of the sample, as the kernel is numerically less than full rank.
  f = t(gaussSamp(Sigma=K, numSamps=1))

  dev.new()
  a = plot(f, xlab='n', ylab= bquote(f[n]))
  dev.copy2eps(file = paste(path, filename, '.eps', sep='')) ## Save plot as eps
  ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
  system(paste('eps2png ', path,filename,'.eps', sep=''))

  return (list(K=K, f=f, x=x))
}