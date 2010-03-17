gpdisimCreate <- function(Ngenes, Ntf, times, y, yvar, options, annotation=NULL, genes = NULL) {

  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")

  if ( Ngenes != dim(y)[2] - 1 )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")

  y <- as.matrix(y)
##   if (options$includeNoise)
##     yvar <- 0 * y
##   else
  yvar <- as.matrix(yvar)

  model <- list(type="gpdisim", y=as.array(y), yvar=as.array(yvar))
  model$includeNoise <- options$includeNoise
  if (("gaussianInitial" %in% names(options)) && options$gaussianInitial)
    model$gaussianInitial <- TRUE
  else
    model$gaussianInitial <- FALSE

  if ("debug" %in% names(options))
    model$debug <- options$debug

  if ("timeSkew" %in% names(options)) {
    model$timeSkew <- options$timeSkew
    times <- times + model$timeSkew
  }
  else
    model$timeSkew <- 0.0

  model$uniqueT <- sort(unique(times))
  lBounds <- c(min(diff(model$uniqueT)),
               (model$uniqueT[length(model$uniqueT)]-model$uniqueT[1]))
  invWidthBounds <- c(2/(lBounds[2]^2), 2/(lBounds[1]^2))

  kernType1 <- list(type="multi", comp=list())
  tieWidth <- 1
  tieRBFVariance <- 2
  kernType1$comp[[1]] <- list(type="parametric", realType="rbf",
                              options=list(isNormalised=TRUE,
                                inverseWidthBounds=invWidthBounds))
  #kernType1$comp[[1]] <- "rbf"
  for ( i in seq(1, Ngenes) ) {
    if (model$gaussianInitial)
      kernType1$comp[[i+1]] <- list(type="parametric", realType="disim",
                                    options=list(gaussianInitial=TRUE,
                                      isNormalised=TRUE,
                                      inverseWidthBounds=invWidthBounds))
    else
      kernType1$comp[[i+1]] <- list(type="parametric", realType="disim",
                                    options=list(isNormalised=TRUE,
                                      inverseWidthBounds=invWidthBounds))
  }
  tieParam <- list(tieDelta="di_decay", tieWidth="inverseWidth",
                   tieSigma="di_variance", tieRBFVariance="rbf.?_variance")

  if ("fixedBlocks" %in% names(options))
    kernType1 <- list(type="parametric", realType=kernType1,
                      options=list(fixedBlocks=options$fixedBlocks))

  if (model$includeNoise) {
    kernType2 <- list(type="multi", comp=array())
    for ( i in seq(1, Ngenes+1) )
      kernType2$comp[i] <- "white"        

    if ("singleNoise" %in% names(options) && options$singleNoise) {
      tieParam$tieNoise <- "white._variance"
    }

    model$kern <- kernCreate(times,
                             list(type="cmpnd", comp=list(kernType1, kernType2)))
    simMultiKernName <- 'model$kern$comp[[1]]'
    simMultiKern <- model$kern$comp[[1]]
  } else {
    model$kern <- kernCreate(times, kernType1)
    simMultiKernName <- 'model$kern'
    simMultiKern <- model$kern
  }
  
  model$kern <- modelTieParam(model$kern, tieParam)

  if (model$includeNoise) {
    for ( i in seq(1, Ngenes+1) ) {
      model$kern$comp[[2]]$comp[[i]]$variance <- 1e-2
    }
  }
  
  model$delta <- 10
  model$sigma <- 1
  for ( i in seq(2, simMultiKern$numBlocks) ) {
    eval(parse(text=paste(simMultiKernName, '$comp[[i]]$di_decay <- model$delta', sep="")))
    eval(parse(text=paste(simMultiKernName, '$comp[[i]]$di_variance <- model$sigma^2', sep="")))
    model$D[i-1] <- simMultiKern$comp[[i]]$decay
    model$S[i-1] <- sqrt(simMultiKern$comp[[i]]$variance)
  }

  set.seed(1)

  yArray <- array(y[,2:(Ngenes+1)], dim = c(dim(y)[1], Ngenes))

  model$numParams <- Ngenes + model$kern$nParams
  model$numGenes <- Ngenes
  model$mu <- apply(yArray, 2, mean)
  model$B <- model$D*yArray[1, ]
  model$m <- array(model$y, length(model$y))
  model$t <- times
  model$realt <- times - model$timeSkew

  model$optimiser <- options$optimiser

  if ( "bTransform" %in% names(options) ) {
    model$bTransform <- options$bTransform
  } else {
    model$bTransform <- "positive"
  }

  if ( !is.null(annotation) ) {
    model$annotation <- annotation
  }

  if (!is.null(genes)) {
    model$genes <- genes
  }

  if ( "fix" %in% names(options) ) {
    params <- gpdisimExtractParam(model, only.values=FALSE)
    model$fix <- options$fix
    if (! "index" %in% names(model$fix)) {
      for ( i in seq(along=model$fix$names) ) {
        J <- grep(model$fix$names[i], names(params))
        if (length(J) != 1) {
          stop("gpdisimCreate: inconsistent fixed parameter specification")
        }
        model$fix$index[i] <- J
      }
    }
  }

  params <- gpdisimExtractParam(model)
  model <- gpdisimExpandParam(model, params)

  return (model)
}


gpdisimDisplay <- function(model, spaceNum=0)  {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Gaussian process driving input single input motif model:\n")
  cat(spacing)
  cat(c("  Number of time points: ", dim(model$t)[1], "\n"), sep="")
  cat(spacing)
  cat(c("  Driving TF: ", model$genes[[1]], "\n"), sep="")
  cat(spacing)
  cat(c("  Target genes (", model$numGenes, "):\n"), sep="")
  for (i in seq(2, length(model$genes))) {
    cat(spacing)
    cat(c("    ", model$genes[[i]], "\n"), sep="")
  }

  if(any(model$mu!=0)) {
    cat(spacing)
    cat("  Basal transcription rate:\n")
    for(i in seq(along=model$mu)) {
      cat(spacing);
      cat(c("    Gene ", i, ": ", model$B[i], "\n"), sep="")
    }
  }
  cat(spacing)
  cat("  Kernel:\n")
  kernDisplay(model$kern, 4+spaceNum)
}
  

gpdisimExtractParam <- function (model, only.values=TRUE) {
  return (gpsimExtractParam(model, only.values=only.values))
}



gpdisimExpandParam <- function (model, params) {

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

  funcName <- optimiDefaultConstraint(model$bTransform)
  # Note: ignores funcName$hasArgs
  func <- get(funcName$func, mode="function")

  model$B <- func(params[(endVal+1):length(params)], "atox")

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern
  
  model$delta <- simMultiKern$comp[[2]]$di_decay
  model$sigma <- simMultiKern$comp[[2]]$di_variance

  for ( i in seq(2, simMultiKern$numBlocks) ) {
    model$D[i-1] <- simMultiKern$comp[[i]]$decay
    model$S[i-1] <- sqrt(simMultiKern$comp[[i]]$variance)
  }

  model$mu <- model$B/model$D
  model <- .gpsimUpdateKernels(model)

  ind <- seq(along=model$t)
  lengthObs <- length(ind)
  model$m[ind] <- model$y[ind]
  ind <- ind + lengthObs
  for ( i in seq(length=model$numGenes) ) {
    model$m[ind] <- model$y[ind]-model$mu[i]*array(1, lengthObs, 1)
    ind <- ind+lengthObs
  }

  return (model)
}




gpdisimObjective <- function (params, model) {
  model <- gpdisimExpandParam(model, params)
  f <- -gpdisimLogLikelihood(model)
  return (f)
}



gpdisimLogLikelihood <- function (model) {
  dim <- length(model$y)

  ll <- -dim*log(2*pi) - model$logDetK - t(model$m) %*% model$invK %*% model$m
  ll <- 0.5*ll

  ## prior contributions
  #if ( "bprior" %in% names(model) ) {
  #  ll <- ll + kernPriorLogProb(model$kern)
  #  ll <- ll + priorLogProb(model$bprior, model$B)
  #}
  return (ll)
}



gpdisimGradient <- function (params, model, ...) {
  model <- gpdisimExpandParam(model, params)
  g <- - gpdisimLogLikeGradients (model)

  return (g)
}



gpdisimLogLikeGradients <- function (model) {

  covGrad <- -model$invK + model$invK %*% model$m %*% t(model$m) %*% model$invK
  covGrad <- 0.5*covGrad

  if ( "proteinPrior" %in% names(model) ) {
    g <- kernGradient(model$kern, model$timesCell, covGrad)
  } else {
    g <- kernGradient(model$kern, model$t, covGrad)
  }

  #if ( "bprior" %in% names(model) ) {
  #  g <- g + kernPriorGradient(model$kern)
  #}

  gmuFull <- t(model$m) %*% model$invK

  if ( "proteinPrior" %in% names(model) ) {
    if ( model$includeNoise ) {
      ind <- model$kern$comp[[1]]$diagBlockDim[1]+(1:model$kern$comp[[1]]$diagBlockDim[2])
      gmu <- array(0, model$numGenes)

      for ( i in seq(length=model$numGenes) ) {
        gmu[i] <- sum(gmuFull[ind])
        ind <- ind + model$kern$comp[[1]]$diagBlockDim[i+1]
      }
    } else {
      ind <- model$kern$diagBlockDim[1]+(1:model$kern$diagBlockDim[2])
      gmu <- array(0, model$numGenes)

      for ( i in seq(length=model$numGenes) ) {
        gmu[i] <- sum(gmuFull[ind])
        ind <- ind + model$kern$diagBlockDim[i+1]
      }
    }
  } else {
    numData <- length(model$t)
    ind <- 1:numData
    ind <- ind + numData
    gmu <- array(0, model$numGenes)
    for ( i in seq(length=model$numGenes) ) {
      gmu[i] <- sum(gmuFull[ind])
      ind <- ind + numData
    }
  }

  gb <- gmu/model$D

  funcName <- optimiDefaultConstraint(model$bTransform)
  # Note: ignores funcName$hasArgs
  func <- get(funcName$func, mode="function")

  ## prior contribution
  #if ( "bprior" %in% names(model) ) {
  #  gb <- gb + priorGradient(model$bprior, model$B)
  #}

  gb <- gb*func(model$B, "gradfact")

  ## gradient for D
  gd <- -gmu*model$B/(model$D*model$D)
  decayIndices <- 5

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern

  for ( i in seq(3, simMultiKern$numBlocks, length=(simMultiKern$numBlocks-2)) )
    decayIndices <- c(decayIndices, tail(decayIndices, 1)+2)

  g[decayIndices] <- g[decayIndices]+gd*expTransform(model$D, "gradfact")

  g <- c(g, gb)

  if ( "fix" %in% names(model) )
    g[model$fix$index] <- 0

  return (g)
}



cgpdisimLogLikeGradients <- function (model) {
  g <- gpdisimLogLikeGradients(model$comp[[1]])

  for ( i in seq(2, length(model$comp), length=(length(model$comp)-1)) )
    g <- g + gpdisimLogLikeGradients(model$comp[[i]])

  return (g)
}



cgpdisimLogLikelihood <- function (model) {
  ll <- 0
  for ( i in seq(along=model$comp) )
    ll <- ll + gpdisimLogLikelihood(model$comp[[i]])

  return (ll)
}



cgpdisimExtractParam <- function (model, only.values=TRUE) {
  return (gpdisimExtractParam(model$comp[[1]], only.values=only.values))
}



cgpdisimExpandParam <- function (model, params) {
  for ( i in seq(along=model$comp) )
    model$comp[[i]] <- gpdisimExpandParam(model$comp[[i]], params)

  return (model)
}



cgpdisimObjective <- function (params, model, trace=0) {
  model <- cgpdisimExpandParam(model, params)
  ll <- - cgpdisimLogLikelihood (model)

  if ( trace & !is.na(ll) )
    cat("Negative Log-Likelihood: ", ll, "\n")

  return (ll)
}



cgpdisimGradient <- function (params, model, ...) {
  model <- cgpdisimExpandParam(model, params)
  g <- - cgpdisimLogLikeGradients (model)

  return (g)
}



cgpdisimUpdateProcesses <- function (model, predt=NULL) {
  for ( i in seq(along=model$comp) )
    model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]], predt=predt)

  return (model)
}



gpdisimUpdateProcesses <- function (model, predt=NULL) {

  if ( "proteinPrior" %in% names(model) ) {
    t <- model$timesCell[[2]]
  } else {
    t <- model$t
  }

  numData <- length(t)
  if (is.null(predt)) {
    predt <- c(seq(min(t), max(t), length=100), t)
  }
  else {
    predt <- c(predt + model$timeSkew, t)
  }

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern

  if (model$gaussianInitial) {
    proteinKern <- kernCreate(model$t, list(type="parametric", realType="sim",
                                            options=list(gaussianInitial=TRUE,
                                              isNormalised=TRUE)))
    proteinKern$initialVariance <- simMultiKern$comp[[2]]$initialVariance
  }
  else
    proteinKern <- kernCreate(model$t, list(type="parametric", realType="sim",
                                            options=list(isNormalised=TRUE)))
  #proteinKern <- kernCreate(model$t, "sim")
  proteinKern$inverseWidth <- simMultiKern$comp[[1]]$inverseWidth
  proteinKern$decay <- model$delta
  proteinKern$variance <- simMultiKern$comp[[2]]$di_variance
  inputKern <- kernCreate(model$t, list(type="parametric", realType="rbf",
                                        options=list(isNormalised=TRUE)))
  inputKern$inverseWidth <- simMultiKern$comp[[1]]$inverseWidth
  inputKern$variance <- 1
  K <- simXrbfKernCompute(proteinKern, inputKern, predt, model$t)
  K <- K * simMultiKern$comp[[1]]$variance
  for ( i in seq(2, length=(simMultiKern$numBlocks-1)) )
    K <- cbind(K,t(disimXsimKernCompute(simMultiKern$comp[[i]], proteinKern, model$t, predt)))

  ymean <- t(matrix(c(0, model$mu), model$numGenes+1, numData))

  predF <- K %*% model$invK %*% c(model$y-ymean)
  #varF <- kernDiagCompute(proteinKern, predt) - apply(t(K)*(model$invK %*% t(K)), 2, sum)
  varF <- simMultiKern$comp[[1]]$variance * kernDiagCompute(proteinKern, predt) - apply(t(K)*(model$invK %*% t(K)), 2, sum) + 1e-15

  Kxx <- multiKernCompute(simMultiKern, predt, model$t)
  meanPredX <- t(matrix(c(0, model$B/model$D), model$numGenes+1, length(predt)))
  predX <- c(meanPredX) + Re(Kxx %*% model$invK %*% c(model$y-ymean))
  varX <- Re( kernDiagCompute(simMultiKern, predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) ) + 1e-13

  predExprs <- matrix(predX, length(predt), model$numGenes+1)
  meanExprs <- t(matrix(apply(predExprs, 2, mean), model$numGenes+1, length(predt)))
  scaleExprs <- t(matrix(sd(predExprs)/sd(as.matrix(model$y)), model$numGenes+1, length(predt)))
  predExprs <- predExprs[1:(length(predt)-length(t)),]

  varExprs <- matrix(varX, length(predt), model$numGenes+1)
  varExprs <- varExprs/scaleExprs/scaleExprs
  varExprs <- varExprs[1:(length(predt)-length(t)),]

  predF <- predF[1:(length(predt)-length(t))]
  varF <- varF[1:(length(predt)-length(t))]
  predt <- predt[1:(length(predt)-length(t))]

  model$predt <- predt - model$timeSkew
  model$predF <- predF
  #model$varF <- abs(varF)
  model$varF <- varF

  model$ypred <- predExprs
  #model$ypredVar <- abs(varExprs)
  model$ypredVar <- varExprs

  return (model)
}
