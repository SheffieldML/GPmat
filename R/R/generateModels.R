generateModels <- function(preprocData = NULL, y = NULL, yvar = NULL, allGenes = NULL, times = NULL, scoreList = NULL, params = NULL, genes = NULL, useGPsim = NULL) {

  options(error = recover)
  options <- list(includeNoise=0, optimiser="CG")
  options$fix$index <- 4
  options$fix$value <- expTransform(c(1, 1), "xtoa")

  if (is.scoreList(scoreList)) {
    params <- scoreList@params
    genes <- scoreList@genes
    useGPsim <- scoreList@useGPsim
  }

  if (!is.null(preprocData)) {
    y <- preprocData$y
    yvar <- preprocData$yvar
    allGenes <- preprocData$genes
    times <- preprocData$times
  }

  listOfModels <- list()
  Nrep <- length(y)

  # for each model in scoreList
  for (i in 1:length(genes)) {

    currentGenes <- genes[[i]]
    currentParams <- params[[i]]
    useGPsimCurrently <- useGPsim[[i]]

    data <- searchExpressionData(y, yvar, allGenes, times, currentGenes)
    currentY <- data$y
    currentYvar <- data$yvar
    currentGenes <- data$genes

    model <- list()

    if (useGPsimCurrently) {
      for (j in 1:Nrep) {
        comp <- gpsimCreate(Ngenes = length(currentGenes), Ntf = 1, times, currentY[[j]], currentYvar[[j]], options)
        model$comp[[j]] <- gpsimExpandParam(comp, currentParams)
      }
    }
    else {
      for (j in 1:Nrep) {
        comp <- gpdisimCreate(Ngenes = length(currentGenes) - 1, Ntf = 1, times, currentY[[j]], currentYvar[[j]], options)
        model$comp[[j]] <- gpdisimExpandParam(comp, currentParams)
      }
    }

    listOfModels[[i]] <- model
  }

  return (listOfModels)
}


searchExpressionData <- function(y, yvar, allGenes, times, searchedGenes) {

  amountOfRows <- length(times)
  Nrep <- length(y)

  # counting the amount of found genes

  # counter for found genes
  k <- 0

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(allGenes)) {
      if (!is.na(charmatch(searchedGene, allGenes[j]))) {
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

  # resetting the counter
  k <- 0

  # copying the data related to the found genes

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(allGenes)) {
      if (!is.na(charmatch(searchedGene, allGenes[j]))) {
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

  newData <- list(y = foundY, yvar = foundYvar, genes = foundGenes, times = times)
  return (newData)
}
