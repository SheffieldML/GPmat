rm(list = ls(all=TRUE))
source("gpsimLoadData.R")
source("miscFunctions.R")

source("kernFunctions.R")
source("mlpKernFunctions.R")
source("optimiFunctions.R")
source("gpsimMapFunctions.R")

options(error = recover)

expType <- "demElkMap"
expNo <- 3

trainingData <- gpsimLoadElkData(3)
y <- trainingData$y
yvar <- trainingData$yvar
genes <- trainingData$genes
times <- trainingData$times
scale <- trainingData$scale
rm(trainingData)

Nrep <- length(y)
Ngenes <- length(genes)
Ntf <- 1

options <- list(includeNoise=1, optimiser="CG")
options$proteinPrior <- list(values=array(-5), times=array(0))
options$kern <- "mlp"
options$nonLinearity <- "activation"
options$startPoint <- 0
options$endPoint <- 8
options$times <- times
options$intPoints <- (options$endPoint-options$startPoint)/0.1+1
options$gParam <- matrix(runif(Ngenes), 1, Ngenes)  ## nrow=ngParam, ncol=numGenes

options$fix$index <- c(3)
options$fix$value <- expTransform(c(1), "xtoa")

model <- list(type="cgpsimMap")
for ( i in seq(length=Nrep) ) {
  options$S <- array(1, Ngenes)
  set.seed(1)
  options$D <- runif(Ngenes)
  mu <- apply(y[[i]], 2, mean)
  options$B <- options$D*mu

  display <- TRUE
  sopt <- TRUE
  options <- gpsimMapInitParam(y[[i]], yvar[[i]], options, display, sopt)
  repNames <- names(model$comp)
  model$comp[[i]] <- gpsimMapCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
  names(model$comp) <- c(repNames, paste("rep", i, sep=""))
  if ( options$kern == "mlp" ) {
    model$comp[[i]]$kern$weightVariance <- 30
    model$comp[[i]]$kern$biasVariance <- 20
    params <- gpsimMapExtractParam(model$comp[[i]])
    model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], params)
  }
}

param <- 0
Nrep <- length(model$comp)
paramvec <- list()
browser()
for ( i in seq(length=Nrep) ) {
  option <- optimiFdefaultOptions()
  model$comp[[i]] <- gpsimMapUpdateF(model$comp[[i]], option)
  fmean <- mean(model$comp[[i]]$f)
  ypredMean <- apply(model$comp[[i]]$ypred, 2, mean)
  ypredScale <- sqrt(apply(model$comp[[i]]$ypred[model$comp[[i]]$timesIndex,], 2, var))
  model$comp[[i]]$B <- (model$comp[[i]]$B - model$comp[[i]]$D*(ypredMean - ypredScale*apply(model$comp[[i]]$y, 2, mean)))/ypredScale
  
  for ( j in seq(length=Ngenes) ) {
    if ( model$comp[[i]]$B[j] < 0 ) {
      model$comp[[i]]$B[j] <- 1e-6
    }
  }
  
  model$comp[[i]]$S <- model$comp[[i]]$S/ypredScale 
  f <- gpsimMapFunctionalExtractParam(model$comp[[i]])
  model$comp[[i]] <- gpsimMapFunctionalExpandParam(model$comp[[i]], f)
  
  paramvec[[i]] <- gpsimMapExtractParam(model$comp[[i]]) 
  param <- param + paramvec[[i]]    

  paramvec <- gpsimMapExtractParam(model$comp[[i]], 2)
  paramVals <- paramvec$values
  param <- param+paramVals
}

param <- param/Nrep

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
save(model, expType, expNo, genes, scale, file=fileName)

displayOption <- 1
gpsimMapElkResults(model, expType, expNo, genes, scale, displayOption)



