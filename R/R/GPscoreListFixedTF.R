GPscoreListFixedTF <- function(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL) {

  # The variable useGPsim determines whether GPsim is used in creating the 
  # models. If its value is false, GPdisim is used. GPsim is used if no TF
  # has been specified.
  useGPsim = is.null(TF)

  if (useGPsim && is.null(knownTargets)) error("There are no known targets for GPsim.")

  amountOfKnownTargets <- length(knownTargets)

  genes = c(TF, knownTargets, testTargets)

  # searching for the data of the specified genes
  searchedData <- searchProcessedData(preprocData, genes)

  # creating a model for each combination and counting the log likelihood of
  # each model

  logLikelihoods <- array(dim = length(testTargets) + 1)
  rankedData <- array(list(NULL), length(testTargets) + 1)
  modelParams <- array(list(NA), length(testTargets) + 1)

  GPrankTargets <- knownTargets
  baseLineParameters <- NULL

  if (!is.null(knownTargets)) {
    baseLineData <- formModel(searchedData, TF, GPrankTargets, useGPsim)
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
      returnData <- formModel(searchedData, TF, GPrankTargets, testTargets[i], useGPsim, fixedParams = TRUE, initParams = baseLineParameters, fixComps = 1:5)
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




formModel <- function(preprocData, TF = NULL, GPrankTargets = NULL, testTarget = NULL, useGPsim = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1) {

    if (!is.null(testTarget)) {
      # taking a test target gene
      GPrankTargets[length(GPrankTargets) + 1] <- testTarget
    }

    error1 <- TRUE
    error2 <- TRUE

    tryCatch({
      data <- GPrank(preprocData, TF, GPrankTargets, useGPsim, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
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
  	  data <- GPrank(preprocData, TF, GPrankTargets, useGPsim, randomize = TRUE, fixedParams = fixedParams, initParams = initParams, fixComps = fixComps)
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
