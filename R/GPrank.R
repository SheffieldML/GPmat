GPrank <- function(preprocData, searchedGenes = "100001_at", search = FALSE, randomize = FALSE, addPriors = FALSE) {

  # GPrank forms an optimized model of the desired genes.
  #
  # preprocData: the preprocessed data to be used
  # searchedGenes: the genes that are searched from the preprocessed data if searching is used
  # search: a logical value determining whether the preprocessed data is searched for the data of specific genes

  options(error = recover)

  gpsimLoc <- "."
  pwd <- getwd()
  setwd(gpsimLoc)

  source("./GPDISIM/gpdisimLoadData.R")
  source("./GPDISIM/gpdisimLinearFunctions.R")
  source("./GPDISIM/disimKernFunctions.R")
  source("gpsimLinearFunctions.R")
  source("miscFunctions.R")
  source("optimiFunctions.R")

  source("kernFunctions.R")
  source("rbfKernFunctions.R")
  source("simKernFunctions.R")
  source("whiteKernFunctions.R")
  source("multiKernFunctions.R")
  source("cmpndKernFunctions.R")

  setwd(pwd)

  source("searchProcessedData.R")
  source("gpdisimCreateFixed.R")

  #expType <- "demGpdisimMef2"
  #expNo <- 1

  # The variable search determines whether the specified genes are searched.
  if (search) {
    newData <- searchProcessedData(preprocData, searchedGenes)
    y <- newData$y
    yvar <- newData$yvar
    times <- newData$times
    genes <- newData$genes
    scale <- newData$scale
  }
  else {
    y <- preprocData$y
    yvar <- preprocData$yvar
    times <- preprocData$times
    genes <- preprocData$genes
    scale <- preprocData$scale
  }

  Nrep <- length(y)

  options <- list(includeNoise=0, optimiser="CG")

  options$fix$index <- c(4, 6)
  options$fix$value <- expTransform(c(1, 1), "xtoa")

  if(addPriors) options$addPriors = TRUE

  Ngenes <- length(genes) - 1
  Ntf <- 1

  # initialise the model
  model <- list(type="cgpdisim")
  for ( i in seq(length=Nrep) ) {
    #repNames <- names(model$comp)
    model$comp[[i]] <- gpdisimCreateFixed(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
    #model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
    #names(model$comp) <- c(repNames, paste("rep", i, sep=""))
  }

  if (randomize) {
    a <- modelExtractParam(model)
    I <- a==0
    n <- length(a)
    a <- array(rnorm(n^2), dim = c(1, n))
    a[I] <- 0
    model <- modelExpandParam(model, a)
  }

  optOptions <- optimiDefaultOptions()

  optOptions$maxit <- 300
  optOptions$optimiser <- "SCG"

  paramvec <- list()
  param <- 0
  for ( i in seq(length=Nrep) ) {
    paramvec[[i]] <- gpdisimExtractParam(model$comp[[i]]) 
    param <- param + paramvec[[i]]    	
  }

  param <- param/Nrep

  optimResult <- SCGoptim(param, fn=cgpdisimObjective, grad=cgpdisimGradient, optOptions, model)
  MLParams <- optimResult$xmin
  #optimResult <- optim(param, fn=cgpdisimObjective, gr=cgpdisimGradient, model, method="BFGS", control=list(maxit=500,trace=10,REPORT=1,parscale=rep(3e-2, length(param))))
  #MLParams <- optimResult$par

  #fileName <- paste(expType, expNo, ".Rdata", sep="")
  #save(model, expType, expNo, genes, scale, file=fileName)

  # optOptions$maxit <- 3000

  # optOptions$fnscale <- 1e1
  # optOptions$trace <- TRUE

  # model <- modelOptimise(model, optOptions)

  for ( i in seq(length=Nrep) ) {
    model$comp[[i]] <- gpdisimExpandParam(model$comp[[i]], MLParams)
    model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]])
  }

  data <- list(model = model, genes = genes)

  return (data)
}
