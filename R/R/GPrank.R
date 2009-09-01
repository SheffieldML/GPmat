GPrank <- function(preprocData, TF = NULL, targets = NULL, useGPsim = FALSE, randomize = FALSE, addPriors = FALSE, search = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

  # GPrank forms an optimized model of the desired genes.
  #
  # preprocData: the preprocessed data to be used
  # TF: the transcription factor of the model
  # targets: the target genes of the model
  # useGPsim: a logical value determining whether GPsim is used; if false, GPdisim is used
  # search: a logical value determining whether the preprocessed data is searched for the data of specific genes

  options(error = recover)

  # The preprocessed data is searched for the data of the specified genes.

  searchedGenes <- c(TF, targets)
  newData <- searchProcessedData(preprocData, searchedGenes)
  y <- newData$y
  yvar <- newData$yvar
  times <- newData$times
  genes <- newData$genes
  scale <- newData$scale

  Nrep <- length(y)

  options <- list(includeNoise=0, optimiser="CG")

  options$fix$index <- 4
  options$fix$value <- expTransform(c(1, 1), "xtoa")

  if(addPriors) options$addPriors = TRUE

  if (useGPsim) {
    Ngenes <- length(genes)
  }
  else {
    Ngenes <- length(genes) - 1
  }
  Ntf <- 1

  # fixing first output sensitivity to fix the scaling
  if (fixedParams && !is.null(initParams)) {
    I <- which(!is.na(initParams))
    for (k in 1:length(I)) {
      options$fix$index[k+1] <- I[k]
      options$fix$value[k+1] <- initParams[I[k]]
    }
  }

  # initializing the model
  model <- list(type="cgpdisim")
  for ( i in seq(length=Nrep) ) {
    #repNames <- names(model$comp)
    if (useGPsim) {
      model$comp[[i]] <- gpsimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
    }
    else {
      model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
    }
    #model$comp[[i]] <- gpdisimCreateFixed(Ngenes, Ntf, times, y[[i]], yvar[[i]], options)
    #names(model$comp) <- c(repNames, paste("rep", i, sep=""))
    if (fixedParams) {
      model$comp[[i]]$kern <- multiKernFixBlocks(model$comp[[i]]$kern, fixComps)
    }
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
    if (useGPsim) {
      paramvec[[i]] <- gpsimExtractParam(model$comp[[i]])
    }
    else {
      paramvec[[i]] <- gpdisimExtractParam(model$comp[[i]])
    } 
    param <- param + paramvec[[i]]    	
  }

  param <- param/Nrep

  cat (c("\n Optimizing genes", TF, targets, sep=" "))

  if(useGPsim) {
    fn <- cgpsimObjective
    grad=cgpsimGradient
  }
  else {
    fn <- cgpdisimObjective
    grad=cgpdisimGradient
  }

  # optimizing the model
  optimResult <- SCGoptim(param, fn, grad, optOptions, model)
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
    if (useGPsim) {
      model$comp[[i]] <- gpsimExpandParam(model$comp[[i]], MLParams)
      model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]])
    }
    else {
      model$comp[[i]] <- gpdisimExpandParam(model$comp[[i]], MLParams)
      model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]])
    }
  }

  data <- list(model = model, genes = genes)

  return (data)
}
