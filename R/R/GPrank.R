GPrank <- function(preprocData, TF = NULL, targets = NULL, useGPsim = FALSE, randomize = FALSE, addPriors = FALSE, search = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

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



GPrankTargets <- function(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL, filterLimit = 0, useMedians = TRUE) {

  # The variable useGPsim determines whether GPsim is used in creating the 
  # models. If its value is false, GPdisim is used. GPsim is used if no TF
  # has been specified.
  useGPsim = is.null(TF)

  if (useGPsim && is.null(knownTargets)) error("There are no known targets for GPsim.")

  amountOfKnownTargets <- length(knownTargets)

  genes <- c(TF, knownTargets, testTargets)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(preprocData$ratioData, genes, filterLimit, useMedians)

  # searching for the data of the specified genes
  searchedData <- searchProcessedData(preprocData, genes)

  logLikelihoods <- array(dim = length(testTargets) + 1)
  rankedData <- array(list(NULL), length(testTargets) + 1)
  modelParams <- array(list(NA), length(testTargets) + 1)

  baseLineParameters <- NULL

  if (!is.null(knownTargets)) {
    baseLineData <- formModel(searchedData, TF, knownTargets, useGPsim)
    logLikelihoods[1] <- baseLineData$ll
    modelParams[[1]] <- baseLineData$params
    rankedData[[1]] <- baseLineData$data

    parameters <- modelExtractParam(baseLineData$data$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*amountOfKnownTargets+4)] <- parameters[1:(2*amountOfKnownTargets+4)]
    t <- 2 * amountOfKnownTargets + 5
    baseLineParameters[(t+2):(t+1+amountOfKnownTargets)] <- parameters[t:(t+amountOfKnownTargets-1)]
  }

  if (length(testTargets) > 0) {
    for (i in 1:length(testTargets)) {
      returnData <- formModel(searchedData, TF, knownTargets, testTargets[i], useGPsim, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
      if (is.null(knownTargets)) {
        logLikelihoods[i] <- returnData$ll
        modelParams[[i]] <- returnData$params
        rankedData[[i]] <- returnData$data
      }
      else {
        logLikelihoods[i+1] <- returnData$ll
        modelParams[[i+1]] <- returnData$params
        rankedData[[i+1]] <- returnData$data
      }
    }
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  scoreList <- list()
  scoreList$LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedData <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))
  for (i in 1:length(sortedIndices)) {
    sortedData[i] <- rankedData[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }
  scoreList$data <- sortedData
  scoreList$params <- sortedModelParams

  return (scoreList)
}



GPrankTFs <- function(preprocData, TFs = NULL, targets = NULL, filterLimit = 0, useMedians = TRUE) {

  if (is.null(targets)) error("There are no known targets for GPsim.")

  amountOfTargets <- length(targets)

  genes = c(TFs, targets)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(preprocData$ratioData, genes, filterLimit, useMedians)

  # searching for the data of the specified genes
  searchedData <- searchProcessedData(preprocData, genes)

  logLikelihoods <- array(dim = length(TFs) + 1)
  rankedData <- array(list(NULL), length(TFs) + 1)
  modelParams <- array(list(NA), length(TFs) + 1)

  baseLineParameters <- NULL

  baseLineData <- formModel(searchedData, TF = NULL, targets, useGPsim = TRUE)
  logLikelihoods[1] <- baseLineData$ll
  modelParams[[1]] <- baseLineData$params
  rankedData[[1]] <- baseLineData$data

  parameters <- modelExtractParam(baseLineData$data$model)
  baseLineParameters <- array(dim = c(1, length(parameters) + 3))
  baseLineParameters[1:(2*amountOfTargets+4)] <- parameters[1:(2*amountOfTargets+4)]
  t <- 2 * amountOfKnownTargets + 5
  baseLineParameters[(t+2):(t+1+amountOfTargets)] <- parameters[t:(t+amountOfTargets-1)]

  if (length(TFs) > 0) {
    for (i in 1:length(TFs)) {
      returnData <- formModel(searchedData, TF = TFs[i], targets, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
      logLikelihoods[i+1] <- returnData$ll
      modelParams[[i+1]] <- returnData$params
      rankedData[[i+1]] <- returnData$data
    }
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  scoreList <- list()
  scoreList$LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedData <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))
  for (i in 1:length(sortedIndices)) {
    sortedData[i] <- rankedData[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }
  scoreList$data <- sortedData
  scoreList$params <- sortedModelParams

  return (scoreList)
}



formModel <- function(preprocData, TF = NULL, knownTargets = NULL, testTarget = NULL, useGPsim = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

    if (!is.null(testTarget)) {
      # taking a test target gene
      knownTargets[length(knownTargets) + 1] <- testTarget
    }

    error1 <- TRUE
    error2 <- TRUE

    tryCatch({
      data <- GPrank(preprocData, TF, knownTargets, useGPsim, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
      error1 <- FALSE
    }, error = function(ex) {
      cat("Stopped due to an error.\n")
    })

    if (error1) {
      success <- FALSE
      i <- 0
      while (!success && i < 10) {
        tryCatch({
	  cat("Trying again with different parameters.\n")
  	  data <- GPrank(preprocData, TF, knownTargets, useGPsim, randomize = TRUE, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
          success <- TRUE
          error2 <- FALSE
        }, error = function(ex) {
          cat("Stopped due to an error.\n")
        })
        i <- i + 1
      }
    }

    else {
      error2 <- FALSE
    }

    if (error2) {
      logLikelihood <- -Inf
      params <- NA
      rankedData <- NA
    }

    else {
      logLikelihood <- logLikelihood(data$model)
      params <- modelExtractParam(data$model)
      rankedData <- data
    }

  returnData <- list()
  returnData$ll <- logLikelihood
  returnData$data <- rankedData
  returnData$params <- params
  return(returnData)
}



searchProcessedData <- function(preprocData, searchedGenes) {

  y <- preprocData$y
  yvar <- preprocData$yvar
  times <- preprocData$times
  genes <- preprocData$genes
  scale <- preprocData$scale

  amountOfRows <- length(times)
  Nrep <- length(y)

  #counting the amount of found genes

  #counter for found genes
  k <- 0

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(genes)) {
      if (!is.na(charmatch(searchedGene, genes[j]))) {
        k <- k + 1
      }
    }
  }

  foundGenes <- array(dim = k)
  foundY <- list()
  foundYvar <- list()
  for (m in 1:Nrep) {
    foundY[[m]] <- array(dim = c(amountOfRows, k))
    foundYvar[[m]] <- array(dim = c(amountOfRows, k))
  }

  #resetting the counter
  k <- 0

  #copying the data related to the found genes

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(genes)) {
      if (!is.na(charmatch(searchedGene, genes[j]))) {
        k <- k + 1
        foundGenes[k] <- searchedGene
	for (m in 1:Nrep) {
          for (l in 1: amountOfRows) {
            foundY[[m]][l, k] <- y[[m]][l, j]
            foundYvar[[m]][l, k] <- yvar[[m]][l, j]
          }
	}
      }
    }
  }
  newData <- list(y = foundY, yvar = foundYvar, genes = foundGenes, times = times, scale = scale)
  return (newData)
}
