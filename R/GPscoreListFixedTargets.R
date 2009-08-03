GPscoreListFixedTargets <- function(preprocData, targets = "100001_at", searchedGenes = "100001_at", search = FALSE) {

  source("GPDISIM/gpdisimLinearFunctions.R")
  source("logLikelihood.R")
  source("GPrank.R")
  source("searchProcessedData.R")

  # The variable searchOption determines whether the target genes are searched.
  if (search) {
    preprocData <- searchProcessedData(preprocData, searchedGenes)
  }

  genes <- preprocData$genes

  if (length(targets) < 0) stop("Invalid number of genes.")

  GPrankSearchedGenes <- array(dim = length(targets) + 1)

  # setting the targets to be tested
  for (i in 1:length(targets)) {
    GPrankSearchedGenes[i + 1] <- targets[i]
  }

  # list of genes with the targets removed, assuming there are no
  # duplicates of them
  otherGenes <- array(dim = length(genes) - length(targets))

  # counter of location in the list of other genes
  k <- 1

  # truth value for checking whether the gene is one of the targets
  isTarget <- FALSE

  # removing the targets from the list of genes
  for (i in 1:length(genes)) {
    for (j in 1:length(targets)) {
      if (!is.na(charmatch(genes[i], targets[j]))) {
	isTarget <- TRUE
      }
    }
    if (!isTarget) {
      otherGenes[k] <- genes[i]
      k <- k + 1
    }
    isTarget <- FALSE
  }

  # creating a model for each TF and counting the log likelihood of each model

  logLikelihoods <- array(dim = length(otherGenes))
  #rankedData <- array(dim = length(otherGenes))
  #rankedData <- vector(mode = 'list', length = length(otherGenes))
  rankedData <- array(list(NULL), length(otherGenes))
  for (i in 1:length(otherGenes)) {

    # taking the TF
    GPrankSearchedGenes[1] <- otherGenes[i]

    data <- GPrank(preprocData, GPrankSearchedGenes, search = TRUE)
    logLikelihoods[i] <- logLikelihood(data$model)
    rankedData[[i]] <- data
  }

  # Sort the log likelihoods.
  sortedValues <- sort(logLikelihoods, decreasing = TRUE, index.return = TRUE)
  scoreList <- list()
  scoreList$LLs <- sortedValues$x
  sortedIndices <- sortedValues$ix

  # Sort the models based on the log likelihoods.
  sortedData <- array(dim = length(sortedIndices))
  for (i in 1:length(sortedIndices)) {
    sortedData[i] <- rankedData[sortedIndices[i]]
  }
  scoreList$data <- sortedData

  return (scoreList)
}
