GPrank <- function(preprocData, TF = NULL, targets = NULL, useGPsim = FALSE, filterLimit = 0, useMedians = TRUE, randomize = FALSE, addPriors = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

  options(error = recover)

  genes <- c(TF, targets)

  # The preprocessed data is searched for the data of the specified genes.
  newData <- searchProcessedData(preprocData, genes)
  y <- newData$y
  yvar <- newData$yvar
  times <- newData$times
  scale <- newData$scale

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(newData$ratioData, filterLimit, useMedians)

  # The data is searched for the data of the filtered genes.
  newData <- searchProcessedData(newData, genes)
  y <- newData$y
  yvar <- newData$yvar
  times <- newData$times
  scale <- newData$scale

  # testing whether the TF passed filtering
  if (!is.null(TF)) {
    if (!genePassedFiltering(TF, genes)) {
      TF <- NULL
    }
  }

  # counter for the next passed gene
  k <- 1

  passedTargets <- array()

  # testing whether the targets passed filtering
  if (!is.null(targets)) {
    for (i in 1:length(targets)) {
      currentTarget <- targets[i]
      if (genePassedFiltering(currentTarget, genes)) {
        passedTargets[k] <- currentTarget
        k <- k + 1
      }
    }
    targets <- passedTargets
  }

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

  param <- modelExtractParam(model)

  cat (c("\n Optimizing genes", genes, sep=" "))

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

  #return (model)

  data <- list(model = model, genes = genes)
  return (data)
}



GPrankTargets <- function(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL, filterLimit = 1.5, useMedians = TRUE) {

  genes <- c(TF, knownTargets, testTargets)

  if (length(genes) < 2) stop("There are too few genes in input.")

  # searching for the data of the specified genes
  searchedData <- searchProcessedData(preprocData, genes)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(searchedData$ratioData, filterLimit, useMedians)

  if (length(genes) < 2) stop("Too few genes passed the filtering.")

  # searching for the data of the genes that passed the filtering
  searchedData <- searchProcessedData(searchedData, genes)

  # testing whether the TF passed filtering
  if (!is.null(TF)) {
    if (!genePassedFiltering(TF, genes)) {
      TF <- NULL
    }
  }

  # counter for the next passed gene
  k <- 1

  passedKnownTargets <- array()

  # testing whether the known targets passed filtering
  if (!is.null(knownTargets)) {
    for (i in 1:length(knownTargets)) {
      currentTarget <- knownTargets[i]
      if (genePassedFiltering(currentTarget, genes)) {
        passedKnownTargets[k] <- currentTarget
        k <- k + 1
      }
    }
    knownTargets <- passedKnownTargets
  }

  k <- 1

  passedTestTargets <- array()

  # testing whether the test targets passed filtering
  if (!is.null(testTargets)) {
    for (i in 1:length(testTargets)) {
      currentTarget <- testTargets[i]
      if (genePassedFiltering(currentTarget, genes)) {
        passedTestTargets[k] <- currentTarget
        k <- k + 1
      }
    }
    testTargets <- passedTestTargets
  }


  # The variable useGPsim determines whether GPsim is used in creating the 
  # models. If its value is false, GPdisim is used. GPsim is used if no TF
  # has been specified.
  useGPsim = is.null(TF)

  if (useGPsim && is.null(knownTargets)) stop("There are no known targets for GPsim.")

  amountOfKnownTargets <- length(knownTargets)

  logLikelihoods <- array(dim = length(testTargets) + 1)
  rankedData <- array(list(NULL), length(testTargets) + 1)
  modelParams <- array(list(NA), length(testTargets) + 1)

  baseLineParameters <- NULL
  fixedParams <- FALSE

  genes <- list()

  if (!is.null(knownTargets) && length(knownTargets) > 0) {
    GPsimUses <- array(dim = length(testTargets) + 1)
  }
  else {
    GPsimUses <- array(dim = length(testTargets))
  }

  if (!is.null(knownTargets) && length(knownTargets) > 0) {
    baseLineData <- formModel(searchedData, TF, knownTargets, useGPsim)
    logLikelihoods[1] <- baseLineData$ll
    modelParams[[1]] <- baseLineData$params
    rankedData[[1]] <- baseLineData$data
    genes[[1]] <- c(TF, knownTargets)
    GPsimUses[1] <- useGPsim

    parameters <- modelExtractParam(baseLineData$data$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*amountOfKnownTargets+4)] <- parameters[1:(2*amountOfKnownTargets+4)]
    t <- 2 * amountOfKnownTargets + 5
    baseLineParameters[(t+2):(t+1+amountOfKnownTargets)] <- parameters[t:(t+amountOfKnownTargets-1)]

    fixedParams <- TRUE
  }

  if (length(testTargets) > 0) {
    for (i in 1:length(testTargets)) {
      returnData <- formModel(searchedData, TF, knownTargets, testTargets[i], useGPsim, fixedParams, initParams = baseLineParameters, fixComps = 1:5)
      if (is.null(knownTargets)) {
        logLikelihoods[i] <- returnData$ll
        modelParams[[i]] <- returnData$params
        rankedData[[i]] <- returnData$data
        genes[[i]] <- c(TF, knownTargets, testTargets[i])
        GPsimUses[i] <- useGPsim
      }
      else {
        logLikelihoods[i+1] <- returnData$ll
        modelParams[[i+1]] <- returnData$params
        rankedData[[i+1]] <- returnData$data
        genes[[i+1]] <- c(TF, knownTargets, testTargets[i])
        GPsimUses[i+1] <- useGPsim
      }
    }
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedData <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))

  for (i in 1:length(sortedIndices)) {
    sortedData[i] <- rankedData[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }

  list <- list()
  list$data <- sortedData

  list$scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)

  return (list)
}



GPrankTFs <- function(preprocData, TFs = NULL, targets = NULL, filterLimit = 1.5, useMedians = TRUE) {

  if (is.null(targets)) stop("There are no targets for GPsim.")

  amountOfTargets <- length(targets)

  genes = c(TFs, targets)

  # searching for the data of the specified genes
  searchedData <- searchProcessedData(preprocData, genes)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(searchedData$ratioData, filterLimit, useMedians)

  if (length(genes) < 2) stop("Too few genes passed the filtering.")

  # searching for the data of the genes that passed the filtering
  searchedData <- searchProcessedData(searchedData, genes)

  # counter for the next passed gene
  k <- 1

  passedTFs <- array()

  # testing whether the TFs passed filtering
  if (!is.null(TFs)) {
    for (i in 1:length(TFs)) {
      currentTF <- TFs[i]
      if (genePassedFiltering(currentTF, genes)) {
        passedTFs[k] <- currentTF
        k <- k + 1
      }
    }
    TFs <- passedTFs
  }

  if (is.null(TFs)) stop("There are no transcription factors for GPsim.")

  k <- 1

  passedTargets <- array()

  # testing whether the targets passed filtering
  if (!is.null(targets)) {
    for (i in 1:length(targets)) {
      currentTarget <- targets[i]
      if (genePassedFiltering(currentTarget, genes)) {
        passedTargets[k] <- currentTarget
        k <- k + 1
      }
    }
    targets <- passedTargets
  }

  if (is.null(targets)) stop("There are no targets for GPsim.")

  logLikelihoods <- array(dim = length(TFs) + 1)
  rankedData <- array(list(NULL), length(TFs) + 1)
  modelParams <- array(list(NA), length(TFs) + 1)

  baseLineParameters <- NULL
  initParams <- FALSE

  amountOfTargets <- length(targets)

  genes <- list()

  if (amountOfTargets > 1) {
    GPsimUses <- array(dim = length(TFs) + 1)
  }
  else {
    GPsimUses <- array(dim = length(TFs))
  }

  if (amountOfTargets > 1) {

    baseLineData <- formModel(searchedData, TF = NULL, targets, useGPsim = TRUE)
    logLikelihoods[1] <- baseLineData$ll
    modelParams[[1]] <- baseLineData$params
    rankedData[[1]] <- baseLineData$data
    genes[[1]] <- targets
    GPsimUses[[1]] <- TRUE

    parameters <- modelExtractParam(baseLineData$data$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*amountOfTargets+4)] <- parameters[1:(2*amountOfTargets+4)]
    t <- 2 * amountOfTargets + 5
    baseLineParameters[(t+2):(t+1+amountOfTargets)] <- parameters[t:(t+amountOfTargets-1)]

    initParams <- TRUE
  }

  for (i in 1:length(TFs)) {
    returnData <- formModel(searchedData, TF = TFs[i], targets, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
    logLikelihoods[i+1] <- returnData$ll
    modelParams[[i+1]] <- returnData$params
    rankedData[[i+1]] <- returnData$data
    genes[[i+1]] <- c(TFs[i], targets)
    GPsimUses[[i+1]] <- FALSE
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedData <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))
  for (i in 1:length(sortedIndices)) {
    sortedData[i] <- rankedData[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }

  list <- list()
  list$data <- sortedData

  list$scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)

  return (list)
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



genePassedFiltering <- function(testedGene, approvedGenes) {

  for (i in 1:length(approvedGenes)) {
    if (!is.na(charmatch(testedGene, approvedGenes[i]))) {
      return (TRUE)
    }
  }

  return (FALSE)
}



searchProcessedData <- function(preprocData, searchedGenes) {

  y <- preprocData$y
  yvar <- preprocData$yvar
  times <- preprocData$times
  genes <- preprocData$genes
  scale <- preprocData$scale
  ratioData <- preprocData$ratioData

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
  foundRatioData <- list()
  foundRatioData$ratios <- list()
  foundRatioData$means <- list()
  foundRatioData$medians <- list()

  for (m in 1:Nrep) {
    foundY[[m]] <- array(dim = c(amountOfRows, k))
    foundYvar[[m]] <- array(dim = c(amountOfRows, k))
    foundRatioData$ratios[[m]] <- array(dim = c(amountOfRows, k))
    foundRatioData$means[[m]] <- array(dim = c(k))
    foundRatioData$medians[[m]] <- array(dim = c(k))
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
	  foundRatioData$means[[m]][k] <- ratioData$means[[m]][j]
	  foundRatioData$medians[[m]][k] <- ratioData$medians[[m]][j]
          for (l in 1: amountOfRows) {
            foundY[[m]][l, k] <- y[[m]][l, j]
            foundYvar[[m]][l, k] <- yvar[[m]][l, j]
            foundRatioData$ratios[[m]][l, k] <- ratioData$ratios[[m]][l, j]
          }
	}
      }
    }
  }

  foundRatioData$genes <- foundGenes

  newData <- list(y = foundY, yvar = foundYvar, genes = foundGenes, times = times, scale = scale, ratioData = foundRatioData)
  return (newData)
}
