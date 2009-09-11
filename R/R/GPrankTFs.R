GPrankTFs <- function(preprocData, TFs = NULL, targets = NULL) {

  if (is.null(targets)) error("There are no known targets for GPsim.")

  amountOfTargets <- length(targets)

  genes = c(TFs, targets)

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
