## DEMGPSAMPLE Simple demonstration of sampling from a covariance function.
demGpSample <- function( bw=FALSE ) {
  require(fields)

  x = as.matrix(seq(-1,1,length=25))
  kern = kernCreate(x, 'rbf')
  kern$inverseWidth = 10

  graphics.off() ## kill all devices
  par(pty="s") #figure(1)
  par(mar=c(4, 4, 2, 6)) ## increase right margin for legend
  K = kernCompute(kern, x)
  colormap=heat.colors(256)
  if (bw)
    colormap=gray.colors(256)

  image(1:dim(K)[2], 1:dim(K)[1] , K[dim(K)[2]:1, ], xlab = "n", ylab = "n")
  image.plot(legend.only = TRUE, nlevel = 256, zlim = c(min(K),max(K)), col = colormap) ## Colorbar
  dev.copy2eps(file = "~/mlprojects/gp/R/inst/doc/gpCovariance.eps") ## Save plot as eps
  ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
  system('eps2png ~/mlprojects/gp/R/inst/doc/gpCovariance.eps')

#   if (exists('printDiagram') && exists('var') && exists(printDiagram))
#     printPlot('gpCovariance', '../tex/diagrams/', '../html') ## !!!

  ## Need to take the real part of the sample, as the kernel is numerically less than full rank.
  f = t(Re(gaussSamp(K, 1)))

  dev.new() #figure(2)
  a = plot(f, xlab='n', ylab= bquote(f[n]))
  dev.copy2eps(file = "~/mlprojects/gp/R/inst/doc/gpSample.eps") ## Save plot as eps
  ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
  system('eps2png ~/mlprojects/gp/R/inst/doc/gpSample.eps')

#   if (exists('printDiagram') && exists('var') && exists(printDiagram))
#     printPlot('gpSample', '../inst/doc/', '../html') ## !!!

  # save(K, x, f, file='demGpSample.RData')
  return (list(K=K, f=f, x=x))
}