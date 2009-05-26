rm(list = ls(all=TRUE))
source("gpsimLoadData.R")
source("miscFunctions.R")

source("kernFunctions.R")
source("mlpKernFunctions.R")
source("optimiFunctions.R")
source("gpsimMapFunctions.R")
source("translateKernFunctions.R")
source("cmpndKernFunctions.R")

expType <- "demEcoliMap"
expNo <- 1

options(error = recover)

trainingData <- gpsimLoadEcoliFullData()
y <- trainingData$y
genes <- trainingData$genes
times <- trainingData$times
scale <- trainingData$scale
rm(trainingData)

Nrep <- length(y)
Ngenes <- length(genes)
Ntf <- 1

yvar <- list()
for ( i in seq(length=Nrep) ) {
  yvar[[i]] <- matrix(0, dim(y[[i]])[1], dim(y[[i]])[2])
}

numGenes <- length(genes)


options <- list(includeNoise=1, optimiser="SCG")
# options$proteinPrior <- list(values=array(), times=array())
options$kern <- c("translate", "mlp")
options$nonLinearity <- "repression"
options$startPoint <- 0
options$endPoint <- 60
options$intPoints <- (options$endPoint-options$startPoint)/0.5+1
options$gParam <- matrix(1, 1, Ngenes) 
options$times <- times
options$bTransform <- "exp"
options$alphaTransform <- ""

options$S <- array(1, Ngenes)
options$D <- runif(Ngenes)
mu <- apply(y[[i]], 2, mean)
options$B <- options$D*mu

options0 <- options

model <- list(type="cgpsimMap")
for ( i in seq(length=Nrep) ) {
  display <- FALSE
  options <- gpsimMapInitParam(y[[i]], yvar[[i]], options0, display)
  if ( options$kern[2] == "mlp" ) {
    options$fix$index <- c(3)
    options$fix$value <- expTransform(c(1), "xtoa")
  } else if ( options$kern[2] == "rbf" ) {
    options$fix$index <- c(2, 19, 20, 21)
    options$fix$value <- expTransform(c(1, 1.75, 0.23, 0.01), "xtoa")
  }

  repNames <- names(model$comp)
  model$comp[[i]] <- gpsimMapCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
  names(model$comp) <- c(repNames, paste("rep", i, sep=""))
  if ( options$kern[2] == "mlp" ) {
    model$comp[[i]]$kern$comp[[1]]$weightVariance <- 30
    model$comp[[i]]$kern$comp[[1]]$biasVariance <- 20
    model$comp[[i]]$kern$centre <- -15
    params <- gpsimMapExtractParam(model$comp[[i]])
    model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], params)
  }
}

paramvec <- list()
param <- 0
for ( i in seq(length=Nrep) ) {
  options <- optimiFdefaultOptions()
  model$comp[[i]] <- gpsimMapUpdateF(model$comp[[i]], options)
  fmean <- mean(model$comp[[i]]$f)
  ypredMean <- apply(model$comp[[i]]$ypred, 2, mean)
  ypredScale <- sqrt(apply(model$comp[[i]]$ypred[model$comp[[i]]$timesIndex,], 2, var))
  model$comp[[i]]$B <- (model$comp[[i]]$B - model$comp[[i]]$D*(ypredMean - ypredScale*mean(model$comp[[i]]$y)))/ypredScale
    
  if ( any(grep("alpha", names(model))) ) 
    model$comp[[i]]$alpha <- model$comp[[i]]$alpha/ypredScale
   
  for ( j in seq(length=numGenes) ) {
    if ( model$comp[[i]]$B[j] < 0 ) {
      model$comp[[i]]$alpha[j] <- model$comp[[i]]$alpha[j] + model$comp[[i]]$B[j]
      model$comp[[i]]$B[j] <- 1e-6
    }
  }
  
  model$comp[[i]]$S <- model$comp[[i]]$S/ypredScale 
  f <- gpsimMapFunctionalExtractParam(model$comp[[i]])
  model$comp[[i]] <- gpsimMapFunctionalExpandParam(model$comp[[i]], f)
  
  paramvec[[i]] <- gpsimMapExtractParam(model$comp[[i]]) 
  param <- param + paramvec[[i]]    
}

param <- param/Nrep;
ll0 <- gpsimMapObjective(param, model)
iters <- 300;

optOptions <- optimiDefaultOptions()
optOptions$maxit <- 300
optOptions$optimiser <- "SCG"

optimResult <- SCGoptim(param, fn=gpsimMapObjective, grad=gpsimMapGradients, optOptions, model)

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, expType, expNo, genes, scale, file=fileName)

for ( i in seq(length=Nrep) ) {
  optimiFoptions <- optimiFdefaultOptions()
  model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], optimResult$xmin)
  model$comp[[i]] <- gpsimMapUpdateF(model$comp[[i]], optimiFoptions)
  model$comp[[i]] <- gpsimMapUpdateYpredVar(model$comp[[i]])
}

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, expType, expNo, genes, scale, file=fileName)

displayOption <- 1
gpsimMapResults(model, expType, expNo, genes, scale, displayOption)




