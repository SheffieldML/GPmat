rm(list = ls(all=TRUE))
source("gpsimLoadData.R")
source("miscFunctions.R")

source("kernFunctions.R")
source("mlpKernFunctions.R")
source("optimiFunctions.R")
source("gpsimMapFunctions.R")

expType <- "demBarencoMap"
expNo <- 1

trainingData <- gpsimLoadBarencoPUMAData()
y <- trainingData$y
yvar <- trainingData$yvar
genes <- trainingData$genes
times <- trainingData$times
scale <- trainingData$scale
rm(trainingData)

Nrep <- length(y)
Ngenes <- length(genes[[1]])
Ntf <- 1

options <- list(includeNoise=1, optimiser="CG")
options$proteinPrior <- list(values=array(-5), times=array(0))
options$kern <- "mlp"
options$nonLinearity <- "activation"
options$startPoint <- 0
options$endPoint <- 12
options$intPoints <- (options$endPoint-options$startPoint)/0.5+1
options$gParam <- matrix(c(1.0129, 0.9605, 0.8504, 0.5950, 0.9337), 1, Ngenes) ## matrix(runif(Ngenes), 1, Ngenes)  ## nrow=ngParam, ncol=numGenes

options$fix$index <- c(18, 3)
options$fix$value <- expTransform(c(0.8, 1), "xtoa")

model <- list(type="cgpsimMap")
for ( i in seq(length=Nrep) ) {
  options$S <- array(1, Ngenes)
  options$D <- c(0.7092, 0.1160, 0.0781, 0.3693, 0.0336)   ## runif(Ngenes)
  mu <- apply(y[[i]], 2, mean)
  options$B <- options$D*mu
  
  repNames <- names(model$comp)
  model$comp[[i]] <- gpsimMapCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
  names(model$comp) <- c(repNames, paste("rep", i, sep=""))
  if ( options$kern == "mlp" ) {
    model$comp[[i]]$kern$weightVariance <- 30
    model$comp[[i]]$kern$biasVariance <- 1000
    params <- gpsimMapExtractParam(model$comp[[i]])
    model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], params)
  }
}

param <- 0
Nrep <- length(model$comp)

for ( i in 1 ) { ## seq(length=Nrep) ) {
  paramvec <- gpsimMapExtractParam(model$comp[[i]], 2)
  paramVals <- paramvec$values
  param <- param+paramVals
}

## param <- param/Nrep

optOptions <- optimiDefaultOptions()
optOptions$maxit <- 300
optOptions$optimiser <- "SCG"

optimResult <- SCGoptim(param, fn=gpsimMapObjective, grad=gpsimMapGradients, optOptions, model)

for ( i in seq(length=Nrep) ) {
  optimiFoptions <- optimiFdefaultOptions()
  model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], optimResult$xmin)
  model$comp[[i]] <- gpsimMapUpdateF(model$comp[[i]], optimiFoptions)
  model$comp[[i]] <- gpsimMapUpdateYpredVar(model$comp[[i]])
}

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, expType, expNo, scale, file=fileName)

displayOption <- 1
gpsimMapBarencoResults(model, expType, expNo, scale, displayOption)



