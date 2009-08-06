gpdisimCreateFixed <- function(Ngenes, Ntf, times, y, yvar, options, annotation=NULL) {

  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")

  if ( Ngenes != dim(y)[2] - 1 )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")

  y <- as.matrix(y)
  yvar <- as.matrix(yvar)
  model <- list(type="gpdisim", y=as.array(y), yvar=as.array(yvar))

  kernType1 <- list(type="multi", comp=array())
  kernType2 <- list(type="multi", comp=array())
  tieWidth <- 1
  tieRBFVariance <- 2
  kernType1$comp[1] <- "rbf"
  for ( i in 1:Ngenes ) {
    kernType1$comp[i+1] <- "disim"
    if ( i==1 ) {
      tieDelta = 3
      tieWidth = c(tieWidth, 4)
      tieSigma = 5
      tieRBFVariance = c(tieRBFVariance, 8)
    } else {
      tieDelta = c(tieDelta, tieDelta[length(tieDelta)]+6)
      tieWidth = c(tieWidth, tieWidth[length(tieWidth)]+6)
      tieSigma = c(tieSigma, tieSigma[length(tieSigma)]+6)
      tieRBFVariance = c(tieRBFVariance, tieRBFVariance[length(tieRBFVariance)]+6)
    }
  }

  tieParam <- list(tieDelta=tieDelta, tieWidth=tieWidth, tieSigma=tieSigma, tieRBFVariance=tieRBFVariance)

  model$kern <- kernCreate(times, kernType1)
  model$kern <- modelTieParam(model$kern, tieParam)

  model$delta <- 10
  model$sigma <- 1

  for ( i in seq(2, model$kern$numBlocks) ) {
    model$kern$comp[[i]]$di_decay <- model$delta
    model$kern$comp[[i]]$di_variance <- model$sigma^2
    model$D[i-1] <- model$kern$comp[[i]]$decay
    model$S[i-1] <- sqrt(model$kern$comp[[i]]$variance)
  }

  set.seed(1)

  yArray <- array(y[,2:(Ngenes+1)], dim = c(dim(y)[1], Ngenes))

  model$numParams <- Ngenes + model$kern$nParams
  model$numGenes <- Ngenes
  model$mu <- apply(yArray, 2, mean)
  model$B <- model$D*yArray[1, ]
  model$m <- array(model$y, length(model$y))
  model$t <- times

  model$optimiser <- options$optimiser

  if ( any(grep("fix", names(options))) ) {
    model$fix <- options$fix
  }

  if ( any(grep("bTransform",names(options))) ) {
    model$bTransform <- options$bTransform
  } else {
    model$bTransform <- "positive"
  }

  if ( !is.null(annotation) ) {
    model$annotation <- annotation
  }

  params <- gpdisimExtractParam(model)
  model <- gpdisimExpandParam(model, params)

  return (model)

}
