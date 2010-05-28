options(error = recover)

library(gpsim)

expType <- "demGpdisimMef2"
expNo <- 1

trainingData <- gpdisimLoadMef2Data()

y <- trainingData$y
yvar <- trainingData$yvar
genes <- trainingData$genes
times <- trainingData$times
scale <- trainingData$scale
rm(trainingData)

genenames <- genes

options <- list(includeNoise=0, optimiser="CG")

options$fix$index <- c(2, 6)
options$fix$value <- expTransform(c(1, 1), "xtoa")

Nrep <- length(y)
Ngenes <- dim(y[[1]])[2] - 1
Ntf <- 1

model <- list(type="cgpdisim")
for ( i in seq(length=Nrep) ) {
  repNames <- names(model$comp)
  model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
  names(model$comp) <- c(repNames, paste("rep", i, sep=""))
}


optOptions <- optimiDefaultOptions()
optOptions$maxit <- 3000

optOptions$fnscale <- 1e1
optOptions$trace <- TRUE

model <- modelOptimise(model, optOptions)

fileName <- paste(expType, expNo, ".Rdata", sep="")
save(model, expType, expNo, genes, scale, file=fileName)

for ( i in seq(length=Nrep) )
  model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]])

gpdisimMef2Results(model, scale, expType, expNo)
