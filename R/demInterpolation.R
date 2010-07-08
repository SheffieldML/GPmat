# DEMINTERPOLATION Demonstrate Gaussian processes for interpolation.
# FORMAT
# DESC runs a simple one-D Gaussian process displaying errorbars.
#
# SEEALSO : gpCreate, demRegression
# 
# COPYRIGHT : Neil D. Lawrence, 2006, 2008
# modifications: Alfredo Kalaitzis, 2010

source('~/mlprojects/gp/R/gpLoadToolkit.R')
gpLoadToolkit()

options = gpOptions(); options$kern$comp = list('rbf','white')

## Create data set
l=9; x = seq(-1, 1, length=l)
trueKern = kernCreate(x, 'rbf')
K = kernCompute(trueKern, x)
## Sample some true function values.
yTrue = gaussSamp(K, 1)

## Create a test set
steps = 2^c(round(log2(l-1)):0); s=0
indTrain=list(); length(indTrain)=length(steps)
indTrain = lapply(indTrain, function(x) x=(seq(1,l,by=steps[(s<<-s+1)]))) # <<-transcends local scope
n = ceiling(sqrt(2*length(steps)-1)); layout.show(layout(matrix(c(1:n^2), n, n, byrow=T)))

figNo = 1
for (i in 1:length(indTrain)) {
  yTrain = as.matrix(yTrue[indTrain[[i]]])
  xTrain = as.matrix(x[indTrain[[i]]])
  kern = kernCreate(x, 'rbf')
  kern$inverseWidth = 5 ## Change inverse variance (1/(lengthScale^2)))

  xTest = as.matrix(seq(-2, 2, length=200))

  Kx = kernCompute(kern, xTest, xTrain)
  Ktrain = kernCompute(kern, xTrain, xTrain)

  invKtrain = .jitCholInv(Ktrain, silent=TRUE)$invM
  yPred = Kx%*%invKtrain%*%yTrain
  yVar = kernDiagCompute(kern, xTest) - rowSums(Kx%*%invKtrain * Kx)

  model = gpCreate(dim(xTrain)[2], dim(yTrain)[2], xTrain, yTrain, options)
  gpPlot(model, xTest, yPred, yVar, ylim=c(-3,3), col='black')

  if (exists('printDiagram') && exists('printDiagram')) 
    printPlot(paste('demInterpolation',as.character(figNo),sep=''),'../tex/diagrams', '../html')
  figNo = figNo + 1

  if (i < length(indTrain)) {
    gpPlot(model, xTest, yPred, yVar, ylim=c(-3,3), col='black')
    diffs = setdiff(indTrain[[i+1]], indTrain[[i]])
    newx = x[diffs]; newy = yTrue[diffs]
    points(newx, newy, pch = 4, cex = 1.5, lwd=3, col='blue')

    if (exists('printDiagram') && exists('printDiagram')) 
      printPlot(paste('demInterpolation',as.character(figNo),sep=''),'../tex/diagrams', '../html')
    figNo = figNo + 1
  }
}
