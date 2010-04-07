GPLearn <- function(preprocData, TF = NULL, targets = NULL,
                    useGpdisim = FALSE, randomize = FALSE,
                    addPriors = FALSE, fixedParams = FALSE,
                    initParams = NULL, initialZero = TRUE,
                    fixComps = NULL, dontOptimise = FALSE,
                    allowNegativeSensitivities = FALSE, quiet = FALSE,
                    gpsimOptions = NULL, allArgs = NULL) {

  if (!is.null(allArgs)) {
    for (i in seq(along=allArgs))
      assign(names(allArgs)[[i]], allArgs[[i]])
  }

  if (is.list(targets))
    targets <- unlist(targets)
  
  if (useGpdisim)
    genes <- c(TF, targets)
  else
    genes <- targets

  # The preprocessed data is searched for the data of the specified genes.
  newData <- .getProcessedData(preprocData[genes,])
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

  if (!initialZero && useGpdisim)
    options$timeSkew <- 1000.0
  #options$gaussianInitial <- TRUE

  if (useGpdisim) {
    Ngenes <- length(genes) - 1
    options$fix$names <- "di_variance"
    options$fix$value <- expTransform(c(1), "xtoa")
  }
  else {
    Ngenes <- length(genes)
    if (allowNegativeSensitivities) {
      options$fix$names <- "sim1_sensitivity"
      options$fix$value <- 1
    }
    else {
      options$fix$names <- "sim1_variance"
      options$fix$value <- expTransform(c(1), "xtoa")
    }
  }
  Ntf <- 1

  if (initialZero && !useGpdisim) {
    options$proteinPrior <- list(values=array(0), times=array(0))

    ## Set the variance of the latent function to 1.
    options$fix$names <- append(options$fix$names, "rbf1_variance")
    options$fix$value <- append(options$fix$value, expTransform(1, "xtoa"))
  }
  
  # fixing first output sensitivity to fix the scaling
  if (fixedParams && !is.null(initParams)) {
    if (is.list(initParams))
      initParams <- unlist(initParams)
    if (!is.null(names(initParams))) {
      options$fix$names <- names(initParams)
      options$fix$value <- initParams
    }
    else {
      I <- which(!is.na(initParams))
      for (k in 1:length(I)) {
        options$fix$index[k+1] <- I[k]
        options$fix$value[k+1] <- initParams[I[k]]
      }
    }      
  }

  if (! is.null(fixComps))
    options$fixedBlocks <- fixComps

  if (allowNegativeSensitivities)
    options$isNegativeS <- TRUE
  
  # initializing the model
  if (useGpdisim)
    model <- list(type="cgpdisim")
  else
    model <- list(type="cgpsim")

  for ( i in seq(length=Nrep) ) {
    #repNames <- names(model$comp)
    if (useGpdisim) {
      model$comp[[i]] <- gpdisimCreate(Ngenes, Ntf, times, t(y[[i]]), t(yvar[[i]]), options, genes = genes)
    }
    else {
      model$comp[[i]] <- gpsimCreate(Ngenes, Ntf, times, t(y[[i]]), t(yvar[[i]]), options, genes = genes)
    }
    #names(model$comp) <- c(repNames, paste("rep", i, sep=""))
    #if (fixedParams) {
    #  model$comp[[i]]$kern <- multiKernFixBlocks(model$comp[[i]]$kern, fixComps)
    #}
  }

  if (randomize) {
    a <- modelExtractParam(model)
    #I <- a==0
    n <- length(a)
    a <- array(rnorm(n), dim = c(1, n))
    #a[I] <- 0
    model <- modelExpandParam(model, a)
  }

  if (!is.null(initParams)) {
    if (!is.list(initParams) && all(is.finite(initParams))
        && is.null(names(initParams)))
      model <- modelExpandParam(model, initParams)
    else if (!is.null(names(initParams))) {
      params <- modelExtractParam(model, only.values=FALSE)
      for (i in seq(along=initParams)) {
        j <- grep(names(initParams)[i], names(params))
        if (length(j) > 0)
          params[j] <- initParams[i]
        else
          warning(paste("Ignoring invalid initial parameter specification:", initParams$names[i]))
      }
      model <- modelExpandParam(model, params)
    }
  }

  optOptions <- optimiDefaultOptions()

  optOptions$maxit <- 3000
  optOptions$optimiser <- "SCG"
  if (quiet)
    optOptions$display <- FALSE

  if (!dontOptimise) {
    message(c("Optimising genes", genes, "\n"), sep=" ")
    model <- modelOptimise(model, optOptions)
  }

  model <- modelUpdateProcesses(model)

  model$allArgs <- list(TF=TF, targets=targets, useGpdisim=useGpdisim,
                        randomize=randomize, addPriors=addPriors,
                        fixedParams=fixedParams, initParams=initParams,
                        initialZero=initialZero, fixComps=fixComps,
                        dontOptimise=dontOptimise, gpsimOptions=gpsimOptions)

  return (new("GPModel", model))
}



GPRankTargets <- function(preprocData, TF = NULL, knownTargets = NULL,
                          testTargets = NULL, filterLimit = 1.8,
                          returnModels = FALSE, options = NULL,
                          scoreSaveFile = NULL) {

  if (is.null(testTargets))
    testTargets <- featureNames(preprocData)

  if (is.list(testTargets))
    testTargets <- unlist(testTargets)
  
  if (is.list(knownTargets))
    knownTargets <- unlist(knownTargets)
  
  if ('var.exprs' %in% assayDataElementNames(preprocData))
    testTargets <- .filterTargets(preprocData[testTargets,], filterLimit)
  else
    testTargets <- featureNames(preprocData[testTargets,])

  if (length(testTargets) < 1)
    stop("No test targets passed filtering")
  
  # The variable useGpdisim determines whether GPDISIM is used in creating the 
  # models. GPSIM (default) is used if no TF has been specified.
  useGpdisim = !is.null(TF)

  if (!useGpdisim && is.null(knownTargets))
    stop("There are no known targets for GPSIM.")

  if (useGpdisim)
    numberOfKnownGenes <- length(knownTargets) + 1
  else
    numberOfKnownGenes <- length(knownTargets)

  logLikelihoods <- rep(NA, length.out=length(testTargets))
  baselogLikelihoods <- rep(NA, length.out=length(testTargets))
  modelParams <- list()
  modelArgs <- list()
  if (returnModels)
    rankedModels <- list()

  genes <- list()

  if (!is.null(options)) {
    allArgs <- options
    allArgs$useGpdisim <- useGpdisim
  }
  else
    allArgs <- list(useGpdisim=useGpdisim)

  if (!is.null(knownTargets) && length(knownTargets) > 0) {
    baselineModel <- .formModel(preprocData, TF, knownTargets, allArgs=allArgs)
    baselineParameters <- modelExtractParam(baselineModel$model,
                                            only.values=FALSE)
    sharedModel <- list(ll=baselineModel$ll,
                        params=baselineModel$params,
                        args=modelStruct(baselineModel$model)$allArgs)

    allArgs$fixedParams <- TRUE
    allArgs$initParams <- baselineParameters
    allArgs$fixComps <- 1:numberOfKnownGenes
  }
  else {
    baselineModel <- NULL
  }

  for (i in seq(along=testTargets)) {
    returnData <- .formModel(preprocData, TF, knownTargets,
                            testTargets[i], allArgs = allArgs)

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
      modelArgs[[i]] <- modelStruct(returnData$model)$allArgs
      if (returnModels)
        rankedModels[[i]] <- returnData$model
    }
    genes[[i]] <- testTargets[[i]]

    testdata <- preprocData[testTargets[i],]
    testdata$experiments <- rep(1, length(testdata$experiments))
    newData <- .getProcessedData(testdata)
    baselogLikelihoods[i] <- .baselineOptimise(newData$y[[1]], newData$yvar[[1]], list(includeNoise=(any(newData$yvar[[1]]==0))))

    if (!is.null(scoreSaveFile)) {
      scoreList <- new("scoreList", params = modelParams,
                       loglikelihoods = logLikelihoods,
                       baseloglikelihoods = baselogLikelihoods,
                       genes = genes, modelArgs = modelArgs,
                       knownTargets = knownTargets, TF = TF,
                       sharedModel = sharedModel)
      save(scoreList, file=scoreSaveFile)
    }
  }

  scoreList <- new("scoreList", params = modelParams,
                   loglikelihoods = logLikelihoods,
                   baseloglikelihoods = baselogLikelihoods,
                   genes = genes, modelArgs = modelArgs,
                   knownTargets = knownTargets, TF = TF,
                   sharedModel = sharedModel)

  if (returnModels)
    return (list(scores=scoreList, models=rankedModels))
  else
    return (scoreList)
}



GPRankTFs <- function(preprocData, TFs, targets,
                      filterLimit = 1.8, 
                      returnModels = FALSE, options = NULL,
                      scoreSaveFile = NULL) {
  if (is.null(targets)) stop("No targets specified.")

  if (is.list(targets))
    targets <- unlist(targets) 

  numberOfTargets <- length(targets)

  genes = c(TFs, targets)

  # Filtering the genes based on the calculated ratios. If the limit is 0, all genes are accepted.
  if ('var.exprs' %in% assayDataElementNames(preprocData))
    TFs <- .filterTargets(preprocData[TFs,], filterLimit)
  else
    TFs <- featureNames(preprocData[TFs,])

  if (length(TFs) < 1)
    stop("No TFs passed the filtering.")

  logLikelihoods <- rep(NA, length.out=length(TFs))
  modelParams <- list()
  modelArgs <- list()
  if (returnModels)
    rankedModels <- list()

  if (!is.null(options)) {
    allArgs <- options
    allArgs$useGpdisim <- TRUE
  }
  else
    allArgs <- list(useGpdisim=TRUE)

  numberOfTargets <- length(targets)

  genes <- list()

  for (i in 1:length(TFs)) {
    returnData <- .formModel(preprocData, TF = TFs[i], targets, allArgs=allArgs)
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
      modelArgs[[i]] <- modelStruct(returnData$model)$allArgs
      if (returnModels)
        rankedModels[[i]] <- returnData$model
    }
    genes[[i]] <- TFs[[i]]

    if (!is.null(scoreSaveFile)) {
      scoreList <- new("scoreList", params = modelParams,
                       loglikelihoods = logLikelihoods,
                       genes = genes, modelArgs = modelArgs,
                       knownTargets = targets, TF = '(see genes)')
      save(scoreList, file=scoreSaveFile)
    }
  }

  scoreList <- new("scoreList", params = modelParams,
                   loglikelihoods = logLikelihoods,
                   genes = genes, modelArgs = modelArgs,
                   knownTargets = targets, TF = '(see genes)')

  if (returnModels)
    return (list(scores=scoreList, models=rankedModels))
  else
    return (scoreList)
}



.formModel <- function(preprocData, TF = NULL, knownTargets = NULL,
                       testTarget = NULL, allArgs = NULL) {

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
    warning("Stopped due to an error.\n")
  })

  if (error1) {
    success <- FALSE
    i <- 0
    allArgs$randomize <- TRUE
    while (!success && i < 10) {
      tryCatch({
        message("Trying again with different parameters.\n")
        model <- GPLearn(preprocData, TF, targets, allArgs = allArgs)
        success <- TRUE
        error2 <- FALSE
      }, error = function(ex) {
        warning("Stopped due to an error.\n")
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


generateModels <- function(preprocData, scores) {
  models <- list()

  # recreate the models for each gene in the scoreList
  for (i in seq(along=params(scores))) {
    args <- modelArgs(scores)[[i]]
    args$initParams <- params(scores)[[i]]
    args$dontOptimise <- TRUE
    models[[i]] <- GPLearn(preprocData, allArgs=args)
  }

  return (models)
}


.getProcessedData <- function(data) {
  times <- data$modeltime
  experiments <- data$experiments

  y <- assayDataElement(data, 'exprs')
  if ('var.exprs' %in% assayDataElementNames(data))
    yvar <- assayDataElement(data, 'var.exprs')
  else
    yvar <- 0 * y

  scale <- sqrt(rowMeans(y^2))
  scaleMat <- scale %*% array(1, dim = c(1, ncol(data)))
  y <- y / scaleMat
  yvar <- yvar / scaleMat^2

  ylist <- list()
  yvarlist <- list()

  expids <- unique(experiments)
  for (i in seq(along=expids)) {
    ylist[[i]] <- y[,experiments==expids[i]]
    yvarlist[[i]] <- yvar[,experiments==expids[i]]
  }
  times <- times[experiments==expids[1]]

  genes <- featureNames(data)

  newData <- list(y = ylist, yvar = yvarlist, genes = genes, times = times)
  return (newData)
}



.filterTargets <- function(data, filterLimit) {
  y <- assayDataElement(data, 'exprs')
  yvar <- assayDataElement(data, 'var.exprs')

  zScores <- rowMeans(y / sqrt(yvar))
  testTargets <- names(which(zScores > filterLimit))

  return (testTargets)
}
