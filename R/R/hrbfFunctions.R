hrbfCreate <- function(Ngenes, times, y, yvar, options, genes = NULL, annotation = NULL) {

  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")
  
  if ( Ngenes != dim(y)[2] )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")
  
  model <- list(type="hrbf", y=as.array(y), yvar=as.array(yvar))

  if ("includeNoise" %in% names(options))
    model$includeNoise <- options$includeNoise
  else
    model$includeNoise <- FALSE

  model$uniqueT <- sort(unique(times))
  if ("lBounds" %in% names(options))
    lBounds <- options$lBounds
  else
    lBounds <- c(min(diff(model$uniqueT)),
                 (model$uniqueT[length(model$uniqueT)]-model$uniqueT[1]))
  invWidthBounds <- c(1/(lBounds[2]^2), 1/(lBounds[1]^2))

  model$isHierarchical <- TRUE
  model$experimentStructure <- options$experiments
  model$masklen <- dim(as.array(options$experiments))[2]

  if ("isSharedHierVariance" %in% names(options))
    model$isSharedHierVariance <- options$isSharedHierVariance
  else
    model$isSharedHierVariance <- FALSE

  if ("isNegativeS" %in% names(options) && options$isNegativeS)
    isNegativeS <- TRUE
  else
    isNegativeS <- FALSE

  myopts <- list(isNormalised=TRUE, inverseWidthBounds=invWidthBounds,
                 isNegativeS=isNegativeS)

  if ("debug" %in% names(options))
    model$debug <- options$debug

  timesCell <- times

  kernType1 <- .rbfKernelSpec(myopts, exps=model$experimentStructure)

  tieParam <- list(tieWidth="inverseWidth")

  if ("fixedBlocks" %in% names(options)) {
    for (k in seq_along(kernType1$comp))
      kernType1$comp[[k]] <-
        list(type="parametric", realType=kernType1$comp[[k]],
             options=list(fixedBlocks=options$fixedBlocks))
  }

  if (model$includeNoise) {
    kernType2 <- list(type="selproj", comp=list(list(type="white")),
                      options=list(expmask=rep(0, model$masklen)))

    model$kern <- kernCreate(model$masklen+1,
                             list(type="cmpnd", comp=list(kernType1, kernType2)))
  } else {
    model$kern <- kernCreate(model$masklen+1, kernType1)
  }
  
  model$kern <- modelTieParam(model$kern, tieParam)

  if (model$includeNoise) {
    model$kern$comp[[2]]$comp[[1]]$variance <- 1e-2
  }

  ## can be removed
  ## if ( model$includeNoise )
  ##   model$yvar <- array(0, length(yvar))

##   if ( model$includeNoise ) {

##     kernType2 <- list(type="multi", comp=array())
##     for ( i in 1:Ngenes )
##       kernType2$comp[i] <- "white"        

##     if ("singleNoise" %in% names(options) && options$singleNoise) {
##       tieParam$tieNoise <- "white\\d+_variance"
##     }

##     model$kern <- kernCreate(model$masklen+1,
##                              list(type="cmpnd", comp=list(kernType1, kernType2)))
##     simMultiKernName <- 'model$kern$comp[[1]]'
##     simMultiKern <- model$kern$comp[[1]]
##   } else {
##     model$kern <- kernCreate(model$masklen+1, kernType1)
##     simMultiKernName <- 'model$kern'
##     simMultiKern <- model$kern
##   }

  ## This is if we need to place priors on parameters ...
  ## ...
  ## ...

  model$kern <- modelTieParam(model$kern, tieParam)

  model$numParams <- model$kern$nParams

  model$numGenes <- Ngenes

  ##model$m <- c(model$y)
  model$t <- t(t(as.array(times)))

  if ("optimiser" %in% names(options))
    model$optimiser <- options$optimiser
  else
    model$optimiser <- "SCG"

  if (!is.null(genes)) {
    model$genes <- genes
  }

  if (!is.null(annotation)  && annotation != "") {
    model$annotation <- annotation
  }

  if ( "fix" %in% names(options) ) {
    params <- modelExtractParam(model, only.values=FALSE)
    model$fix <- options$fix
    if (! "index" %in% names(model$fix)) {
      for ( i in seq(along=model$fix$names) ) {
        J <- grep(model$fix$names[i], names(params))
        if (length(J) != 1) {
          stop("hrbfCreate: inconsistent fixed parameter specification")
        }
        model$fix$index[i] <- J
      }
    }
  }

  params <- hrbfExtractParam(model)
  model <- hrbfExpandParam(model, params)
  
  return (model)
}


hrbfDisplay <- function(model, spaceNum=0)  {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Gaussian process multiple output model:\n")
  cat(spacing)
  cat(c("  Number of time points: ", dim(model$t)[1], "\n"), sep="")
  cat(spacing)
  cat(c("  Number of genes: ", model$numGenes, "\n"), sep="")
  
  cat(spacing)
  cat("  Kernel:\n")
  kernDisplay(model$kern, 4+spaceNum)
}
  

hrbfExtractParam <- function (model, only.values=TRUE,
                               untransformed.values=FALSE) {
  if (untransformed.values && "fix" %in% names(model)) {
    params <- hrbfExtractParam(model, only.values=TRUE,
                                untransformed.values=FALSE)
    origparams <- params
    for ( i in seq(along=model$fix$index) )
      params[model$fix$index[i]] <- model$fix$value[i]

    if (! all(params==origparams)) {
      model <- hrbfExpandParam(model, params)
    }
  }
  
  params <- kernExtractParam(model$kern, only.values=only.values,
                             untransformed.values=untransformed.values)

  return (Re(params))
}



hrbfExpandParam <- function (model, params) {
  if ( is.list(params) )
    params <- params$values

  params <- Re(params)
  if ( "fix" %in% names(model) ) 
    for ( i in seq(along=model$fix$index) )
      params[model$fix$index[i]] <- model$fix$value[i]

  if ( length(params) != model$numParams )
    stop("Parameter vector is incorrect length.")

  startVal <- 1
  endVal <- model$kern$nParams
  model$kern <- kernExpandParam(model$kern, params[startVal:endVal])

  ##model$mu <- model$B/model$D
  model$mu <- 0

  model <- .hrbfUpdateKernels(model)

##   ind <- seq(along=model$t)
##   lengthObs <- length(ind)
##   for ( i in seq(length=model$numGenes) ) {
##     model$m[ind] <- model$y[ind]-model$mu[i]*array(1, lengthObs, 1)
##     ind <- ind+lengthObs
##   }

  return (model)
}



.hrbfUpdateKernels <- function (model) {
  eps <-  1e-6

  k <- Re(kernCompute(model$kern, cbind(model$experimentStructure, model$t)))
  ##noiseVar <- c(as.array(model$yvar))

  model$K <- k ##+diag(as.array(noiseVar))
  invK <- .jitCholInv(model$K, silent=TRUE)

  if ( is.nan(invK[1]) ) {
    if ("debug" %in% names(model) && model$debug) {
      cat("kern$decay = \n", model$D, "\n")
      cat("\n")
      cat("kern$sensitivity = \n", model$S, "\n")
      cat("\n")
      cat("kern$flength = \n", model$kern$comp[1]$flength, "\n")
      cat("\n")
    }
    stop("Singular chol(K) matrix!")
  }
  
  model$invK <- invK$invM

  if ( invK$jitter > 1e-4 && "debug" %in% names(model) && model$debug )
    warning(paste("Warning: hrbfUpdateKernels added jitter of", signif(invK$jitter, digits=4)))

  model$logDetK <- 2* sum( log ( diag(invK$chol) ) )
 
  return (model)
}



hrbfOptimise <- function (model, options, ...) {
  if ( "optimiser" %in% names(options) ) {
    funcName <- paste(options$optimiser, "optim", sep="")
  } else {
    funcName <- "CGoptim"
  }
  func <- get(funcName, mode="function")
  
  params <- modelExtractParam(model)
  newParams <- func(params, fn=modelObjective, grad=modelGradient, options, model)
  
  model <- hrbfExpandParam(model, newParams$xmin)
  model$llscore <- newParams$objective
  if ( funcName == "CGoptim" )
    model$lnSchFail <- newParams$lnSchFail

  return (model)
}



hrbfObjective <- function (params, model) {
  model <- hrbfExpandParam(model, params)
  f <- -hrbfLogLikelihood(model)
  return (f)
}



hrbfLogLikelihood <- function (model) {
  dim <- dim(model$y)[1]

  ll <- 0
  for (i in seq(dim(model$y)[2]))
    ll <- ll - dim*log(2*pi) - model$logDetK - t(model$y[,i]) %*% model$invK %*% model$y[,i]
  
  ll <- 0.5*ll

  ## prior contributions
  ll <- ll + kernPriorLogProb(model$kern)
  return (ll)
}



hrbfLogLikeGradients <- function (model) {
  covGrad <- array(0, dim=dim(model$invK))
  for (i in seq(dim(model$y)[2]))
    covGrad <- covGrad - model$invK + model$invK %*% model$y[,i] %*% t(model$y[,i]) %*% model$invK
  
  covGrad <- 0.5*covGrad

  g <- kernGradient(model$kern, cbind(model$experimentStructure, model$t), covGrad)
  g <- g + kernPriorGradient(model$kern)

  if ( "fix" %in% names(model) ) 
    g[model$fix$index] <- 0

  return (g)  
}
 

hrbfGradient <- function (params, model, ...) {
  model <- hrbfExpandParam(model, params)
  g <- - hrbfLogLikeGradients (model)

  return (g)
}


hrbfUpdateProcesses <- function (model, predt=NULL) {
  t <- model$t

  return (model)
  ##numData <- length(t)
  if (is.null(predt))
    predt <- c(seq(min(t), max(t), length=100), t)
  ## ymean <- t(matrix(model$y[1,], model$numGenes, numData))

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern

  #meanPredX <- t(matrix(model$B/model$D, model$numGenes, length(predt)))
  meanPredX <- t(matrix(0, model$numGenes, length(predt)))

  Kxx <- multiKernCompute(simMultiKern, predt, t)
  predX <- c(meanPredX) + Re( Kxx%*%model$invK%*%model$m)
  varX <- Re( kernDiagCompute(simMultiKern, predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) ) 

  predExprs <- matrix(predX, length(predt), model$numGenes)
  varExprs <- matrix(varX, length(predt), model$numGenes)
  predExprs <- predExprs[1:(length(predt)-length(t)),]
  varExprs <- varExprs[1:(length(predt)-length(t)),]                

  model$predt <- predt[1:(length(predt)-length(t))]

  model$ypred <- predExprs
  model$ypredVar <- abs(varExprs)

  return (model)
}
