require(Biobase)

GPLearn <- function(preprocData, TF = NULL, targets = NULL,
                    useGpdisim = FALSE, randomize = FALSE,
                    addPriors = FALSE, fixedParams = FALSE,
                    initParams = NULL, initialZero = TRUE,
                    fixComps = NULL, dontOptimise = FALSE,
                    gpsimOptions = NULL, allArgs = NULL) {

  if (!is.null(allArgs)) {
    for (i in seq(along=allArgs))
      assign(names(allArgs)[[i]], allArgs[[i]])
  }
  
  if (useGpdisim)
    genes <- c(TF, targets)
  else
    genes <- targets

  # The preprocessed data is searched for the data of the specified genes.
  newData <- searchProcessedData(preprocData, genes)
  y <- newData$y
  yvar <- newData$yvar
  times <- newData$times

  Nrep <- length(y)

  options <- list(includeNoise=0, optimiser="SCG")

  if (any(yvar[[1]] == 0))
    options$includeNoise = 1
  
  if (!is.null(gpsimOptions)) {
    for (i in seq(along=gpsimOptions))
      options[[names(gpsimOptions)[[i]]]] <- gpsimOptions[[i]]
  }
  
  if (addPriors)
    options$addPriors <- TRUE

  if (!initialZero)
    options$timeSkew <- 1000.0
  #options$gaussianInitial <- TRUE

  if (useGpdisim) {
    Ngenes <- length(genes) - 1
    options$fix$names <- "di_variance"
    options$fix$value <- expTransform(c(1), "xtoa")
  }
  else {
    Ngenes <- length(genes)
    options$fix$names <- "sim1_variance"
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
  if (useGpdisim)
    model <- list(type="cgpdisim")
  else
    model <- list(type="cgpsim")

  for ( i in seq(length=Nrep) ) {
    #repNames <- names(model$comp)
    if (useGpdisim) {
      model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options, genes = genes)
    }
    else {
      model$comp[[i]] <- gpsimCreate(Ngenes, Ntf, times, y[[i]], yvar[[i]], options, genes = genes)
    }
    #names(model$comp) <- c(repNames, paste("rep", i, sep=""))
    if (fixedParams) {
      model$comp[[i]]$kern <- multiKernFixBlocks(model$comp[[i]]$kern, fixComps)
    }
  }

  if (randomize) {
    a <- modelExtractParam(model)
    #I <- a==0
    n <- length(a)
    a <- array(rnorm(n), dim = c(1, n))
    #a[I] <- 0
    model <- modelExpandParam(model, a)
  }

  if (!is.null(initParams) && all(is.finite(initParams)))
    model <- modelExpandParam(model, initParams)

  optOptions <- optimiDefaultOptions()

  optOptions$maxit <- 3000
  optOptions$optimiser <- "SCG"

  cat (c("\n Optimizing genes", genes, sep=" "))

  if (!dontOptimise)
    model <- modelOptimise(model, optOptions)

  model <- modelUpdateProcesses(model)

  model$allArgs <- list(TF=TF, targets=targets, useGpdisim=useGpdisim,
                        randomize=randomize, addPriors=addPriors,
                        fixedParams=fixedParams, initParams=initParams,
                        initialZero=initialZero, fixComps=fixComps,
                        dontOptimise=dontOptimise, gpsimOptions=gpsimOptions)

  return (GPmodel(model))
}



GPrankTargets <- function(preprocData, TF = NULL, knownTargets = NULL,
                          testTargets = NULL, filterLimit = 1.8,
                          returnModels = FALSE, options = NULL,
                          scoreSaveFile = NULL) {

  if (is.null(testTargets))
    testTargets <- preprocData@genes

  if (any(is.na(preprocData@zScores[testTargets])))
    testTargets <- names(preprocData@zScores[testTargets])
  else
    testTargets <- names(which(preprocData@zScores[testTargets] > filterLimit))

  if (length(testTargets) < 1)
    stop("No test targets passed filtering")
  
  # The variable useGpdisim determines whether GPDISIM is used in creating the 
  # models. GPSIM (default) is used if no TF has been specified.
  useGpdisim = !is.null(TF)

  if (!useGpdisim && is.null(knownTargets)) stop("There are no known targets for GPSIM.")

  numberOfKnownTargets <- length(knownTargets)

  logLikelihoods <- rep(NA, length.out=length(testTargets))
  modelParams <- list()
  modelArgs <- list()
  if (returnModels)
    rankedModels <- list()

  baseLineParameters <- NULL
  fixedParams <- FALSE

  genes <- list()

  if (!is.null(options)) {
    allArgs <- options
    allArgs$useGpdisim <- useGpdisim
  }
  else
    allArgs <- list(useGpdisim=useGpdisim)
  
  ## CHECKME!
  if (!is.null(knownTargets) && length(knownTargets) > 0) {
    baseLineData <- formModel(preprocData, TF, knownTargets, allArgs=allArgs)
    baselogLikelihoods[1] <- baseLineData$ll
    basemodelParams[[1]] <- baseLineData$params
    baserankedModels[[1]] <- baseLineData$model
    basegenes[[1]] <- c(TF, knownTargets)

    parameters <- modelExtractParam(baseLineData$model)
    baseLineParameters <- array(dim = c(1, length(parameters) + 3))
    baseLineParameters[1:(2*numberOfKnownTargets+4)] <- parameters[1:(2*numberOfKnownTargets+4)]
    t <- 2 * numberOfKnownTargets + 5
    baseLineParameters[(t+2):(t+1+numberOfKnownTargets)] <- parameters[t:(t+numberOfKnownTargets-1)]

    fixedParams <- TRUE
    allArgs$fixedParams <- TRUE
    allArgs$initParams <- baseLineParameters
    allArgs$fixedComps <- 1:numberOfKnownTargets
  }

  if (length(testTargets) > 0) {
    for (i in 1:length(testTargets)) {
      returnData <- formModel(preprocData, TF, knownTargets, testTargets[i], allArgs = allArgs)

      if (!is.finite(returnData$ll)) {
        logLikelihoods[i] <- NA
        modelParams[[i]] <- NA
        modelArgs[[i]] <- NA
        if (returnModels)
          rankedModels[[i]] <- NA
      }
      else {
        logLikelihoods[i] <- returnData$ll
        modelParams[[i]] <- returnData$params
        modelArgs[[i]] <- returnData$model@model$allArgs
        if (returnModels)
          rankedModels[[i]] <- returnData$model
      }
      genes[[i]] <- testTargets[[i]]

      if (!is.null(scoreSaveFile)) {
        scoreList <- scoreList(params = modelParams, LLs = logLikelihoods, genes = genes, modelArgs = modelArgs, knownTargets = knownTargets, TF = TF)
        save(scoreList, file=scoreSaveFile)
      }
    }
  }

  scoreList <- scoreList(params = modelParams, LLs = logLikelihoods, genes = genes, modelArgs = modelArgs, knownTargets = knownTargets, TF = TF)

  if (returnModels)
    return (list(scores=scoreList, models=rankedModels))
  else
    return (scoreList)
}



GPrankTFs <- function(preprocData, TFs = NULL, targets = NULL, filterLimit = 1.5, returnScoreList = TRUE, returnModels = FALSE) {

  if (is.null(targets)) stop("No targets specified.")

  numberOfTargets <- length(targets)

  genes = c(TFs, targets)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  genes <- filterGenes(preprocData$ratioData, filterLimit)

  if (length(genes) < 2) stop("Too few genes passed the filtering.")

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

    baseLineData <- formModel(preprocData, TF = NULL, targets, useGPsim = TRUE)
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
    returnData <- formModel(preprocData, TF = TFs[i], targets, useGPsim = FALSE, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
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



formModel <- function(preprocData, TF = NULL, knownTargets = NULL, testTarget = NULL, allArgs = NULL) {

  if (!is.null(testTarget))
    targets <- append(knownTargets, testTarget)
  else
    targets <- knownTargets

  error1 <- TRUE
  error2 <- TRUE

  tryCatch({
    model <- GPLearn(preprocData, TF, targets, allArgs = allArgs)
    error1 <- FALSE
  }, error = function(ex) {
    cat("Stopped due to an error.\n")
  })

  if (error1) {
    success <- FALSE
    i <- 0
    allArgs$randomize <- TRUE
    while (!success && i < 10) {
      tryCatch({
        cat("Trying again with different parameters.\n")
        model <- GPLearn(preprocData, TF, targets, allArgs = allArgs)
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
    logLikelihood <- modelLogLikelihood(model)
    params <- modelExtractParam(model)
    rankedModel <- model
  }

  return (list(ll=logLikelihood, model=rankedModel, params=params))
}


generateModels <- function(preprocData, scores, indices=NULL) {
  models <- list()

  if (is.null(indices))
    indices <- seq(along=scores@params)
  
  # recreate the models for each gene in the scoreList
  j <- 1
  for (i in indices) {
    args <- scores@modelArgs[[i]]
    args$initParams <- scores@params[[i]]
    args$dontOptimise <- TRUE
    models[[j]] <- GPLearn(preprocData, allArgs=args)
    j <- j+1
  }

  return (models)
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
  times <- preprocData@times

  y <- list()
  yvar <- list()
  for (i in seq(along=preprocData@y)) {
    y[[i]] <- as.matrix(preprocData@y[[i]][,searchedGenes])
    yvar[[i]] <- as.matrix(preprocData@yvar[[i]][,searchedGenes])
  }

  genes <- preprocData@genes[searchedGenes]

  newData <- list(y = y, yvar = yvar, genes = genes, times = times)
  return (newData)
}
