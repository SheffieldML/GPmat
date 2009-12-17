GPrank <- function(preprocData, TF = NULL, targets = NULL, useGPsim = TRUE, filterLimit = 0, useMedians = TRUE, randomize = FALSE, addPriors = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

  options(error = recover)

  genes <- c(TF, targets)

  # The preprocessed data is searched for the data of the specified genes.
  newData <- searchProcessedData(preprocData, genes)

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

  if(addPriors) options$addPriors = TRUE

  if (useGPsim) {
    Ngenes <- length(genes)
    options$fix$names <- "sim1_variance"
    options$fix$value <- expTransform(c(1), "xtoa")
  }
  else {
    Ngenes <- length(genes) - 1
    options$fix$names <- "di_variance"
    options$fix$value <- expTransform(c(1), "xtoa")
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
  if (useGPsim)
    model <- list(type="cgpsim")
  else
    model <- list(type="cgpdisim")

  for ( i in seq(length=Nrep) ) {
    #repNames <- names(model$comp)
    if (useGPsim) {
      model$comp[[i]] <- gpsimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options, genes = genes)
    }
    else {
      model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options, genes = genes)
    }
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

  model <- modelOptimise(model, optOptions)

  model <- modelUpdateProcesses(model)

  return (model)
}



GPrankTargets <- function(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL, filterLimit = 1.5, useMedians = TRUE, returnScoreList = TRUE, returnModels = FALSE) {

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

  numberOfKnownTargets <- length(knownTargets)

  logLikelihoods <- array(dim = length(testTargets) + 1)
  rankedModels <- array(list(NULL), length(testTargets) + 1)
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
    rankedModels[[1]] <- baseLineData$model
    genes[[1]] <- c(TF, knownTargets)
    GPsimUses[1] <- useGPsim

    parameters <- modelExtractParam(baseLineData$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*numberOfKnownTargets+4)] <- parameters[1:(2*numberOfKnownTargets+4)]
    t <- 2 * numberOfKnownTargets + 5
    baseLineParameters[(t+2):(t+1+numberOfKnownTargets)] <- parameters[t:(t+numberOfKnownTargets-1)]

    fixedParams <- TRUE
  }

  if (length(testTargets) > 0) {
    for (i in 1:length(testTargets)) {
      returnData <- formModel(searchedData, TF, knownTargets, testTargets[i], useGPsim, fixedParams, initParams = baseLineParameters, fixComps = 1:5)
      if (is.null(knownTargets)) {
        logLikelihoods[i] <- returnData$ll
        modelParams[[i]] <- returnData$params
        rankedModels[[i]] <- returnData$model
        genes[[i]] <- c(TF, knownTargets, testTargets[i])
        GPsimUses[i] <- useGPsim
      }
      else {
        logLikelihoods[i+1] <- returnData$ll
        modelParams[[i+1]] <- returnData$params
        rankedModels[[i+1]] <- returnData$model
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
  sortedModels <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))

  for (i in 1:length(sortedIndices)) {
    sortedModels[i] <- rankedModels[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }

  if (returnScoreList && returnModels) {
    data <- list()
    data$scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)
    data$models <- sortedModels
    return (data)
  }

  if (returnScoreList) {
    scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)
    return (scoreList)
  }

  if (returnModels) {
    return (sortedModels)
  }
}



GPrankTFs <- function(preprocData, TFs = NULL, targets = NULL, filterLimit = 1.5, useMedians = TRUE, returnScoreList = TRUE, returnModels = FALSE) {

  if (is.null(targets)) stop("There are no targets for GPsim.")

  numberOfTargets <- length(targets)

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
  rankedModels <- array(list(NULL), length(TFs) + 1)
  modelParams <- array(list(NA), length(TFs) + 1)

  baseLineParameters <- NULL
  initParams <- FALSE

  numberOfTargets <- length(targets)

  genes <- list()

  if (numberOfTargets > 1) {
    GPsimUses <- array(dim = length(TFs) + 1)
  }
  else {
    GPsimUses <- array(dim = length(TFs))
  }

  if (numberOfTargets > 1) {

    baseLineData <- formModel(searchedData, TF = NULL, targets, useGPsim = TRUE)
    logLikelihoods[1] <- baseLineData$ll
    modelParams[[1]] <- baseLineData$params
    rankedModels[[1]] <- baseLineData$model
    genes[[1]] <- targets
    GPsimUses[[1]] <- TRUE

    parameters <- modelExtractParam(baseLineData$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*numberOfTargets+4)] <- parameters[1:(2*numberOfTargets+4)]
    t <- 2 * numberOfTargets + 5
    baseLineParameters[(t+2):(t+1+numberOfTargets)] <- parameters[t:(t+numberOfTargets-1)]

    initParams <- TRUE
  }

  for (i in 1:length(TFs)) {
    returnData <- formModel(searchedData, TF = TFs[i], targets, useGPsim = FALSE, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
    logLikelihoods[i+1] <- returnData$ll
    modelParams[[i+1]] <- returnData$params
    rankedModels[[i+1]] <- returnData$model
    genes[[i+1]] <- c(TFs[i], targets)
    GPsimUses[[i+1]] <- FALSE
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedModels <- array(dim = length(sortedIndices))
  sortedModelParams <- array(dim = length(sortedIndices))
  for (i in 1:length(sortedIndices)) {
    sortedModels[i] <- rankedModels[sortedIndices[i]]
    sortedModelParams[i] <- modelParams[sortedIndices[i]]
  }

  if (returnScoreList && returnModels) {
    data <- list()
    data$scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)
    data$models <- sortedModels
    return (data)
  }

  if (returnScoreList) {
    scoreList <- scoreList(params = sortedModelParams, LLs = LLs, genes = genes, useGPsim = GPsimUses)
    return (scoreList)
  }

  if (returnModels) {
    return (sortedModels)
  }
}



formModel <- function(preprocData, TF = NULL, knownTargets = NULL, testTarget = NULL, useGPsim = TRUE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

    if (!is.null(testTarget)) {
      # taking a test target gene
      knownTargets <- append(knownTargets, testTarget)
    }

    error1 <- TRUE
    error2 <- TRUE

    tryCatch({
      model <- GPrank(preprocData, TF, knownTargets, useGPsim, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
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
  	  model <- GPrank(preprocData, TF, knownTargets, useGPsim, randomize = TRUE, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
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
      rankedModel <- NA
    }

    else {
      logLikelihood <- logLikelihood(model)
      params <- modelExtractParam(model)
      rankedModel <- model
    }

  returnData <- list()
  returnData$ll <- logLikelihood
  returnData$model <- rankedModel
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

  numberOfRows <- length(times)
  Nrep <- length(y)

  #counting the number of found genes

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
    foundY[[m]] <- array(dim = c(numberOfRows, k))
    foundYvar[[m]] <- array(dim = c(numberOfRows, k))
    foundRatioData$ratios[[m]] <- array(dim = c(numberOfRows, k))
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
          for (l in 1: numberOfRows) {
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


setClass("scoreList", 
	representation(params = "list", LLs = "array", genes = "list", useGPsim = "array"))

	#prototype = list(params = list(), LLs = array(), genes = list(), useGPsim = array()))


scoreList <- function(params = list(), LLs = array(), genes = list(), useGPsim = array()) {

  new("scoreList", params = params, LLs = LLs, genes = genes, useGPsim = useGPsim)
}


is.scoreList <- function(object) {
  return (class(object) == "scoreList")
}
