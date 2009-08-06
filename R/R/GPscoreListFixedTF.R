GPscoreListFixedTF <- function(preprocData, TF = "100001_at", searchedGenes = "100001_at", search = FALSE, amountOfTargetGenes = 1) {

  # The variable searchOption determines whether the target genes are searched.
  if (search) {
    preprocData <- searchProcessedData(preprocData, searchedGenes)
  }

  genes <- preprocData$genes

  if (amountOfTargetGenes < 0) stop("Invalid number of genes.")

  GPrankSearchedGenes <- array(dim = c(amountOfTargetGenes + 1))

  # setting the transcription factor to be tested
  GPrankSearchedGenes[1] <- TF

  # list of genes with the transcription factor removed, assuming there are no
  # duplicates of the transcription factor
  otherGenes <- array(dim = length(genes) - 1)

  # counter of location in the list of other genes
  j <- 1

  # removing the transcription factor from the list of genes
  for (i in 1: length(genes)) {
    if (is.na(charmatch(genes[i], TF))) {
      otherGenes[j] <- genes[i]
      j <- j + 1
    }
  }

  # determining possible gene combinations of the remaining genes
  combinations <- combn(otherGenes, amountOfTargetGenes)
  amountOfCombinations <- dim(combinations)[2]

  # creating a model for each combination and counting the log likelihood of
  # each model

  logLikelihoods <- array(dim = amountOfCombinations)
  #rankedData <- array(dim = amountOfCombinations)
  #rankedData <- vector(mode = 'list', length = amountOfCombinations)
  rankedData <- array(list(NULL), amountOfCombinations)
  for (i in 1:amountOfCombinations) {

    # taking the target genes of this combination
    for (j in 1:amountOfTargetGenes) {
      GPrankSearchedGenes[j + 1] <- combinations[j, i]
    }

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
