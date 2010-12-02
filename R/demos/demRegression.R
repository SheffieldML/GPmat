# DEMREGRESSION Demonstrate Gaussian processes for regression.
# FORMAT
# DESC runs a simple one-D Gaussian process displaying errorbars.
#
# SEEALSO : gpCreate, demInterpolation
# 
# COPYRIGHT : Neil D. Lawrence, 2006
# modifications: Alfredo Kalaitzis, 2010
demRegression <- function(filename='demRegression') {
#   source('~/mlprojects/gp/R/gpLoadToolkit.R')
#   gpLoadToolkit()

  noiseLevel = 0.2
  noiseVar = noiseLevel^2
  options = gpOptions(); options$kern$comp = list('rbf','white')

  ## Create data set
  l=9; x = as.matrix(seq(-1, 1, length=l))
  trueKern = kernCreate(x, 'rbf')
  K = kernCompute(trueKern, x) + diag(dim(x)[1])*noiseVar
  ## Sample some true function values.
  yTrue = gaussSamp(Sigma=K, numSamps=1)

  ## Create a test set
  steps = 2^c(round(log2(l-1)):0); s=0
  indTrain=list(); length(indTrain)=length(steps)
  indTrain = lapply(indTrain, function(x) x=(seq(1,l,by=steps[(s<<-s+1)]))) # <<-transcends local scope
#   n = ceiling(sqrt(2*length(steps)-1)); layout.show(layout(matrix(c(1:n^2), n, n, byrow=T)))
  graphics.off()
  dev.new(); plot.new()

  path = '~/mlprojects/gp/R/inst/doc/'
  figNo = 1
  for (i in 1:length(indTrain)) {
    if (i > 0) {
      yTrain = as.matrix(yTrue[indTrain[[i]]])
      xTrain = as.matrix(x[indTrain[[i]]])
      kern = kernCreate(x, 'rbf')
      ## Change inverse variance (1/(lengthScale^2)))
      kern$inverseWidth = 5

      xTest = as.matrix(seq(-2, 2, length=200))

      Kx = kernCompute(kern, xTest, xTrain)
      Ktrain = kernCompute(kern, xTrain, xTrain)

      invKtrain = .jitCholInv(Ktrain + diag(dim(Ktrain)[1])*noiseVar, silent=TRUE)$invM
      yPred = Kx%*%invKtrain%*%yTrain
      yVar = kernDiagCompute(kern, xTest) - rowSums(Kx%*%invKtrain * Kx)

      model = gpCreate(dim(xTrain)[2], dim(yTrain)[2], xTrain, yTrain, options)
      gpPlot(model, xTest, yPred, yVar, ylim=c(-3,3), col='black')
      
      pathfilename = paste(path,filename,'_', figNo, '.eps', sep='')
      dev.copy2eps(file = pathfilename) ## Save plot as eps
      ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
      system(paste('eps2png ', pathfilename, sep=''))

      figNo = figNo + 1
    } else {
      xTest = as.matrix(seq(-2, 2, length=200))
      kern = kernCreate(xTest, 'rbf')
      ## Change inverse variance (1/(lengthScale^2)))
      kern$inverseWidth = 5

      yPred = array(0,dim(xTest))
      yVar = kernDiagCompute(kern, xTest)

      plot(xTest, yPred, col='blue', lwd=.5)
      zeroAxes(xTest)
      
      pathfilename = paste(path,filename,'_', figNo, '.eps', sep='')
      dev.copy2eps(file = pathfilename) ## Save plot as eps
      ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
      system(paste('eps2png ', pathfilename, sep=''))

      figNo = figNo + 1
    }

    if (i < length(indTrain)) {
      gpPlot(model, xTest, yPred, yVar, ylim=c(-3,3), col='black')
      if (i < 1) diffs = indTrain[[i+1]]
      else diffs = setdiff(indTrain[[i+1]], indTrain[[i]])
      newx = x[diffs]; newy = yTrue[diffs]
      points(newx, newy, pch = 4, cex = 1.5, lwd=3, col='blue')

      pathfilename = paste(path,filename,'_', figNo, '.eps', sep='')
      dev.copy2eps(file = pathfilename) ## Save plot as eps
      ## Convert to png. Needs the 'eps2png' facility. If not already installed: 'sudo apt-get install eps2png'
      system(paste('eps2png ', pathfilename, sep=''))

      figNo = figNo + 1
    }
  }

  ## Convert the .png files to one .gif file using ImageMagick. 
  ## The -delay flag sets the time between showing
  ## the frames, i.e. the speed of the animation.
  system(paste('convert -delay 80 ',path,filename,'*.png ', path,filename,'.gif', sep=''))
}