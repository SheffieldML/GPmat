rm(list = ls(all=TRUE))
source("gpsimLoadData.R")
source("gpsimLinearFunctions.R")
source("miscFunctions.R")

source("kernFunctions.R")
source("rbfKernFunctions.R")
source("simKernFunctions.R")
source("whiteKernFunctions.R")
source("multiKernFunctions.R")
source("cmpndKernFunctions.R")
source("optimiFunctions.R")

expType <- "demBarenco"
expNo <- 2

options <- list(includeNoise=1, optimiser="CG")
options$proteinPrior <- list(values=array(0), times=array(0))

if ( any(i <- grep("proteinPrior",names(options))) ) {
  options$fix$index <- c(9, 10, 2)
  options$fix$value <- expTransform(c(0.8, 1, 1), "xtoa")
} else {
  options$fix$index <- c(8, 9)
  options$fix$value <- expTransform(c(0.8, 1), "xtoa")
}

trainingData <- gpsimLoadBarencoMAS5Data()
y <- trainingData$y
yvar <- trainingData$yvar
genes <- trainingData$genes
times <- trainingData$times
scale <- trainingData$scale
rm(trainingData)

Nrep <- length(y)
Ngenes <- length(genes[[1]])
Ntf <- 1

model <- list(type="cgpsim")
for ( i in seq(length=Nrep) ) {
  repNames <- names(model$comp)
  model$comp[[i]] <- gpsimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
  names(model$comp) <- c(repNames, paste("rep", i, sep=""))
}

optOptions <- optimiDefaultOptions()
optOptions$maxit <- 1000
optOptions$fnscale <- 1e1
optOptions$trace <- TRUE

model <- modelOptimise(model, optOptions)

for ( i in seq(length=Nrep) )
  model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]]) 

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, scale, expType, expNo, file=fileName)
gpsimBarencoResults(model, scale, expType, expNo)

