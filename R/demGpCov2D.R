demGpCov2D <- function(ind=c(1,2), bw=FALSE) {
  library(fields)
  graphics.off() ## kill all devices

  source('~/mlprojects/gp/R/gaussSamp.R')
  source('~/mlprojects/gp/R/demGpSample.R')
  K = K[ind, ind]
  f = f[ind]
  x = x[ind]
  print(K)

  dev.new()
  lo = layout(matrix(c(1:4),2,2,byrow=T)); layout.show(lo)

#   dev.new();
  basePlot(K, ind)

#   dev.new() #figure(2)
  basePlot(K, ind)
  cont2 = lines(c(f[1],f[1]), c(-1,1), col='green')

#   dev.new() #figure(3)
  basePlot(K, ind)
  cont2 = lines(c(f[1],f[1]), c(-1,1), col='green')

  ## Compute conditional mean and variance
  f2Mean = K[1, 2]/K[1,1]*f[1]
  f2Var = K[2, 2] - K[1, 2]/K[1, 1]*K[1, 2]
  yval = as.matrix(seq(-1, 1, length=200))
  pdfVal = 1/sqrt(2*pi*f2Var)*exp(-0.5*(yval-f2Mean)*(yval-f2Mean)/f2Var)
  pdf = lines(pdfVal*0.25, yval, col='red')
#   set(pdf, 'linewidth', conditionalSize)
#   set(pdf, 'linestyle', conditionalLineStyle, 'color', conditionalLineColour)

  if (bw )
    app = 'bw'
  else 
    app = ''

  for (figNo in 2:4) {
    dev.set(figNo) #figure(figNo)
    if (exists('printDiagram') && exists('var') && exists(printDiagram))
      printPlot(paste('demGpCov2D',as.character(ind[1]),'_',as.character(ind[2]),'_', 
	as.character(figNo),app,sep=""), '../tex/diagrams', '../html')
  }
}

basePlot <- function (K, ind) {
  # BASEPLOT Plot the contour of the covariance.
  # DESC creates the basic plot.

  eigVecs = eigen(K)
  U = eigVecs$vectors[,dim(eigVecs$vectors)[2]:1]; V = rev(eigVecs$values) ## in incr. order
  V[V<0] = as.complex(V[V<0])
  r = Re(sqrt(V))
  theta = seq(0, 2*pi, length=200)
  xy = cbind(r[1]*sin(theta), r[2]*cos(theta))%*%U
  plot(0, type = "n", xlim=c(min(xy[,1]), max(xy[,1])),
    ylim=c(min(xy[,2]), max(xy[,2])), xlab='', ylab='') ## lines only works on existing plots
  cont = lines(xy[, 1], xy[, 2], col='blue')

#   t = text(0.5, -0.2, ['f_{' num2str(ind(1)) '}']);
#   t =[t text(-0.2, -0.5, ['f_{' num2str(ind(2)) '}'])];

# 
#   set(gca, 'xtick', [-1  0  1])
#   set(gca, 'ytick', [-1 0 1])
# 
#   set(t, 'fontname', 'times')
#   set(t, 'fontsize', 24)
#   set(t, 'fontangle', 'italic')
# 
#   zeroAxes(gca, 0.025, 18, 'times')
  abline(v = 0, h = 0) ## axes and co-centric circles
  int = (max(xy[,1])-min(xy[,1]))/100
  lines(int*cos(seq(0, 2*pi, l=100)), int*sin(seq(0, 2*pi, l=100)))

#   ax = gca

#   return list(ax, cont, t)
}