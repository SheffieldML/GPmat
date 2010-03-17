gpsimCreate <- function(Ngenes, Ntf, times, y, yvar, options, genes = NULL) {

  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")
  
  if ( Ngenes != dim(y)[2] )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")
  
  model <- list(type="gpsim", y=as.array(y), yvar=as.array(yvar))

  kernType1 <- list(type="multi", comp=list())
  tieParam <- list(tieWidth="inverseWidth")

  model$uniqueT <- sort(unique(times))
  lBounds <- c(min(diff(model$uniqueT)),
               (model$uniqueT[length(model$uniqueT)]-model$uniqueT[1]))
  invWidthBounds <- c(2/(lBounds[2]^2), 2/(lBounds[1]^2))

  if ("isNegativeS" %in% names(options) && options$isNegativeS)
    isNegativeS = TRUE
  else
    isNegativeS = FALSE

  if ("debug" %in% names(options))
    model$debug <- options$debug

  ## proteinPrior encodes observation of the latent function.
  if ( "proteinPrior" %in% names(options) ) {
    model$proteinPrior <- options$proteinPrior
    kernType1$comp[[1]] <- list(type="parametric", realType="rbf",
                                options=list(isNormalised=TRUE,
                                  inverseWidthBounds=invWidthBounds))

    if ( "times" %in% names(model$proteinPrior) )
      timesCell <- list(protein=model$proteinPrior$times)
    else
      timesCell <- list(protein=times)

    for ( i in 1:Ngenes ) {
      kernType1$comp[[i+1]] <- list(type="parametric", realType="sim",
                                    options=list(isNormalised=TRUE,
                                      inverseWidthBounds=invWidthBounds,
                                      isNegativeS=isNegativeS))
      timesCell <- c(timesCell, list(mRNA=times))
    }

    model$timesCell <- timesCell
  } else {
    timesCell <- times
    
    kernType1 <- list(type="multi", comp=list())
    for ( i in 1:Ngenes ) {
      kernType1$comp[[i]] <- list(type="parametric", realType="sim",
                                  options=list(isNormalised=TRUE,
                                    inverseWidthBounds=invWidthBounds,
                                    isNegativeS=isNegativeS))
    }
  }

  if ("fixedBlocks" %in% names(options))
    kernType1 <- list(type="parametric", realType=kernType1,
                      options=list(fixedBlocks=options$fixedBlocks))

  model$includeNoise <- options$includeNoise

  ## can be removed
  ## if ( model$includeNoise )
  ##   model$yvar <- array(0, length(yvar))

  if ( model$includeNoise ) {

    kernType2 <- list(type="multi", comp=array())
    if ( "proteinPrior" %in% names(options) ) {
      for ( i in 1:(Ngenes+1) )
        kernType2$comp[i] <- "white"
    } else {
      for ( i in 1:Ngenes )
        kernType2$comp[i] <- "white"        
    }     

    if ("singleNoise" %in% names(options) && options$singleNoise) {
      tieParam$tieNoise <- "white._variance"
    }

    model$kern <- kernCreate(timesCell,
                             list(type="cmpnd", comp=list(kernType1, kernType2)))
    simMultiKernName <- 'model$kern$comp[[1]]'
    simMultiKern <- model$kern$comp[[1]]
  } else {
    model$kern <- kernCreate(timesCell, kernType1)
    simMultiKernName <- 'model$kern'
    simMultiKern <- model$kern
  }

  ## This is if we need to place priors on parameters ...
  ## ...
  ## ...

  model$kern <- modelTieParam(model$kern, tieParam)

  if ( "proteinPrior" %in% names(options) ) {
    for ( i in seq(2,simMultiKern$numBlocks,length.out=simMultiKern$numBlocks-1) ) {
      model$D[i-1] <- simMultiKern$comp[[i]]$decay
      model$S[i-1] <- simMultiKern$comp[[i]]$sensitivity
    }
  } else {
    for ( i in seq(along=simMultiKern$comp) ) {
      model$D[i] <- simMultiKern$comp[[i]]$decay
      model$S[i] <- simMultiKern$comp[[i]]$sensitivity
    }
  }

  model$numParams <- Ngenes + model$kern$nParams
  model$numGenes <- Ngenes
  model$mu <- apply(y, 2, mean)
  model$B <- model$D*model$mu

  if ( "proteinPrior" %in% names(model) ) {
    model$m <- c(model$proteinPrior$values, model$y)
  } else {
    model$m <- c(model$y)
    model$t <- times
  }

  if ( "bTransform" %in% names(options) ) {
    model$bTransform <- options$bTransform
  } else {
    model$bTransform <- "positive"
  }

  model$optimiser <- options$optimiser

  if (!is.null(genes)) {
    model$genes <- genes
  }

  if ( "fix" %in% names(options) ) {
    params <- modelExtractParam(model, only.values=FALSE)
    model$fix <- options$fix
    if (! "index" %in% names(model$fix)) {
      for ( i in seq(along=model$fix$names) ) {
        J <- grep(model$fix$names[i], names(params))
        if (length(J) != 1) {
          stop("gpsimCreate: inconsistent fixed parameter specification")
        }
        model$fix$index[i] <- J
      }
    }
  }

  params <- gpsimExtractParam(model)
  model <- gpsimExpandParam(model, params)
  
  return (model)
}


gpsimDisplay <- function(model, spaceNum=0)  {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Gaussian process single input motif model:\n")
  cat(spacing)
  cat(c("  Number of time points: ", dim(model$t)[1], "\n"), sep="")
  cat(spacing)
  cat(c("  Number of genes: ", model$numGenes, "\n"), sep="")
  
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
  

gpsimExtractParam <- function (model, only.values=TRUE) {
  funcName <- optimiDefaultConstraint(model$bTransform)
  func <- get(funcName$func, mode="function")

  if ( only.values ) {
    params <- kernExtractParam(model$kern)
    # Note: ignores funcName$hasArgs
    params <- c(params, func(model$B, "xtoa"))
  } else {
    params <- kernExtractParam(model$kern, only.values)
    # Note: ignores funcName$hasArgs
    Bparams <- func(model$B, "xtoa")
    for ( i in seq(along=Bparams) ) {
      names(Bparams)[i] <- paste("Basal", i, sep="")
    }
    params <- c(params, Bparams)
  }

  if ( "fix" %in% names(model) ) 
    for ( i in seq(along=model$fix$index) )
      params[model$fix$index[i]] <- model$fix$value[i]

  params <- Re(params)        

  return (params)
}



gpsimExpandParam <- function (model, params) {
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
  func <- get(funcName$func, mode="function")

  # Note: ignores funcName$hasArgs
  model$B <- func(params[(endVal+1):length(params)], "atox")

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern

  if ( "proteinPrior" %in% names(model) ) {
    for ( i in seq(2,model$kern$numBlocks,length.out=model$kern$numBlocks-1) ) {
      model$D[i-1] <- simMultiKern$comp[[i]]$decay
      model$S[i-1] <- simMultiKern$comp[[i]]$sensitivity
    }
  } else {
    for ( i in seq(length.out=model$kern$numBlocks) ) {
      model$D[i] <- simMultiKern$comp[[i]]$decay
      model$S[i] <- simMultiKern$comp[[i]]$sensitivity
    }
  }

  model$mu <- model$B/model$D

  model <- .gpsimUpdateKernels(model)

  if ( "proteinPrior" %in% names(model) ) {
    yInd <- seq(1, simMultiKern$diagBlockDim[2])
    mInd <- simMultiKern$diagBlockDim[1] + yInd

    for ( i in seq(length=model$numGenes) ) {
      model$m[mInd] <- model$y[yInd]-model$mu[i]*array(1, length(yInd), 1)
      yInd <- yInd+simMultiKern$diagBlockDim[i+1]
      mInd <- mInd+simMultiKern$diagBlockDim[i+1]
    }
  } else {
    ind <- seq(along=model$t)
    lengthObs <- length(ind)
    for ( i in seq(length=model$numGenes) ) {
      model$m[ind] <- model$y[ind]-model$mu[i]*array(1, lengthObs, 1)
      ind <- ind+lengthObs
    }
  }

  return (model)
}



.gpsimUpdateKernels <- function (model) {
  eps <-  1e-6

  if ( ("proteinPrior" %in% names(model)) && ("timesCell" %in% names(model)) ) {
    k <- Re(kernCompute(model$kern, model$timesCell))
    noiseVar <- c(array(eps, model$kern$comp[[1]]$diagBlockDim[1], 1), model$yvar)
  } else {
    k <- Re(kernCompute(model$kern, model$t))
    noiseVar <- c(as.array(model$yvar))
  }

  model$K <- k+diag(as.array(noiseVar))
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
    warning(paste("Warning: gpsimUpdateKernels added jitter of", signif(invK$jitter, digits=4)))

  model$logDetK <- 2* sum( log ( diag(invK$chol) ) )
 
  return (model)
}



cgpsimOptimise <- function (model, options, ...) {
  if ( "optimiser" %in% names(options) ) {
    funcName <- paste(options$optimiser, "optim", sep="")
  } else {
    funcName <- "CGoptim"
  }
  func <- get(funcName, mode="function")
  
  params <- modelExtractParam(model)
  newParams <- func(params, fn=modelObjective, grad=modelGradient, options, model)
  
  model <- cgpsimExpandParam(model, newParams$xmin)
  model$llscore <- newParams$objective
  if ( funcName == "CGoptim" )
    model$lnSchFail <- newParams$lnSchFail

##  if ( length(options) < 7 ) {
##    newOptions <- options
##    options <- optimiDefaultOptions()
##    for ( i in seq(along=newOptions) ) 
##      if ( !is.na(newOptions[[i]]) )
##        options[[i]] <- newOptions[[i]]
##  }

##  params <- cgpsimExtractParam(model)

##  options1 <- optimiDefaultOptions()
##  options1$maxit <- 200
##  options1$fnscale <- 1e3
##  options1$trace <- TRUE

##  modelopt <- optim(params, modelObjective, gr=modelGradient, method="CG", control=list(trace=FALSE, type=2, maxit=options1$maxit, fnscale=options1$fnscale, reltol=options1$reltol), hessian=FALSE, model, options1$trace)
  
##  modelopt <- optim(modelopt$par, modelObjective, gr=modelGradient, method="CG", control=list(trace=FALSE, type=2, maxit=options$maxit, fnscale=options$fnscale, reltol=options$reltol), hessian=options$hessian, model, options$trace)

##  model <- cgpsimExpandParam(model, modelopt$par)
  
  return (model)
}



gpsimObjective <- function (params, model) {
  model <- gpsimExpandParam(model, params)
  f <- -gpsimLogLikelihood(model)
  return (f)
}



gpsimLogLikelihood <- function (model) {
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



gpsimLogLikeGradients <- function (model) {
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
    gmu <- array(0, model$numGenes)
    for ( i in seq(length=model$numGenes) ) {
      gmu[i] <- sum(gmuFull[ind])
      ind <- ind + numData
    }    
  }

  gb <- gmu/model$D

  funcName <- optimiDefaultConstraint(model$bTransform)
  func <- get(funcName$func, mode="function")

  ## prior contribution
  #if ( "bprior" %in% names(model) ) {
  #  gb <- gb + priorGradient(model$bprior, model$B)
  #}

  # Note: ignores funcName$hasArgs
  gb <- gb*func(model$B, "gradfact")

  ## gradient for D
  gd <- -gmu*model$B/(model$D*model$D)
  if ( model$kern$numBlocks == 1 ) {
    decayIndices <- 1
  } else if ( model$kern$numBlocks>1 ) {
    if ( "proteinPrior" %in% names(model) ) {
      decayIndices <- 3
      for ( i in seq(3, model$kern$numBlocks, length=(model$kern$numBlocks-2)) )
        decayIndices <- c(decayIndices, decayIndices[length(decayIndices)]+2)
    } else {
      decayIndices <- c(1, 4)
      for ( i in seq(3, model$kern$numBlocks, length=(model$kern$numBlocks-2)) )
        decayIndices <- c(decayIndices, decayIndices[length(decayIndices)]+2)
    }
  }

  g[decayIndices] <- g[decayIndices]+gd*expTransform(model$D, "gradfact")

  g <- c(g, gb)

  if ( "fix" %in% names(model) ) 
    g[model$fix$index] <- 0

  return (g)  
}
 


cgpsimLogLikeGradients <- function (model) {
  g <- gpsimLogLikeGradients(model$comp[[1]])

  for ( i in seq(2, length(model$comp), length=(length(model$comp)-1)) )
    g <- g + gpsimLogLikeGradients(model$comp[[i]])
  
  return (g)
}



cgpsimLogLikelihood <- function (model) {
  ll <- 0
  for ( i in seq(along=model$comp) )
    ll <- ll + gpsimLogLikelihood(model$comp[[i]])

  return (ll)
}



cgpsimExtractParam <- function (model, only.values=TRUE) {
  return (gpsimExtractParam(model$comp[[1]], only.values=only.values))
}



cgpsimExpandParam <- function (model, params) {
  for ( i in seq(along=model$comp) )
    model$comp[[i]] <- gpsimExpandParam(model$comp[[i]], params)

  return (model)
}



cgpsimObjective <- function (params, model, trace=0) {
  model <- cgpsimExpandParam(model, params)
  ll <- - cgpsimLogLikelihood (model)

  if ( trace & !is.na(ll) )
    cat("Negative Log-Likelihood: ", ll, "\n")

  return (ll)
}



cgpsimGradient <- function (params, model, ...) {
  model <- cgpsimExpandParam(model, params)
  g <- - cgpsimLogLikeGradients (model)

  return (g)
}



cgpsimUpdateProcesses <- function (model, predt=NULL) {
  for ( i in seq(along=model$comp) )
    model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]], predt=predt)

  return (model)
}



gpsimGradient <- function (params, model, ...) {
  model <- gpsimExpandParam(model, params)
  g <- - gpsimLogLikeGradients (model)

  return (g)
}


gpsimUpdateProcesses <- function (model, predt=NULL) {
  if ( "proteinPrior" %in% names(model) ) {
    t <- model$timesCell[[2]]
  } else {
    t <- model$t
  }

  ##numData <- length(t)
  if (is.null(predt))
    predt <- c(seq(min(t), max(t), length=100), t)
  ## ymean <- t(matrix(model$y[1,], model$numGenes, numData))

  if (model$includeNoise)
    simMultiKern <- model$kern$comp[[1]]
  else
    simMultiKern <- model$kern

  if ( "proteinPrior" %in% names(model) ) {
    predTimeCell <- list()
    for ( i in seq(length=simMultiKern$numBlocks) )
      predTimeCell[[i]] <- predt

    Kxx <- multiKernCompute(simMultiKern, predTimeCell, model$timesCell)
    diagKxx <- kernDiagCompute(simMultiKern, predTimeCell)
    ## x <- c(model$proteinPrior$values, model$y-ymean)
    ind <- 1:length(predTimeCell[[1]])
    predFull <- list()
    varFull <- list()
    for ( indBlock in seq(length=simMultiKern$numBlocks) ) {
      K <- Kxx[ind,]
      diagK <- diagKxx[ind,]
      predFull[[indBlock]] <- Re( K %*% model$invK %*% model$m )
      varFull[[indBlock]] <- Re( diagK - t(apply(t(K)*(model$invK %*% t(K)), 2, sum)) )
      ind <- ind + length(predTimeCell[[indBlock]])
    }

    meanPredX <- t(matrix(model$B/model$D, model$numGenes, length(predt)))

    predF <- predFull[[1]]
    varF <- varFull[[1]]
    predExprs <- matrix(, nrow=length(predt), ncol=model$numGenes)
    varExprs <- matrix(, nrow=length(predt), ncol=model$numGenes)

    for ( i in seq(length=model$numGenes) ) {
      predExprs[,i] <- meanPredX[,i] + predFull[[i+1]]
      varExprs[,i] <- varFull[[i+1]]
    }

    predExprs <- predExprs[1:(length(predt)-length(t)),]
    varExprs <- varExprs[1:(length(predt)-length(t)),]

  } else {
    proteinKern <- kernCreate(model$t, "rbf")
    proteinKern$inverseWidth <- simMultiKern$comp[[1]]$inverseWidth

    K <- simXrbfKernCompute(simMultiKern$comp[[1]], proteinKern, model$t, predt)
    for ( i in seq(2, length=(simMultiKern$numBlocks-1)) ) 
      K <- rbind(K, simXrbfKernCompute(simMultiKern$comp[[i]], proteinKern, model$t, predt))

    predF <- Re( t(K)%*%model$invK%*%model$m)
    varF <- Re( kernDiagCompute(proteinKern, predt) - apply(K*(model$invK%*%K), 2, sum) )
    meanPredX <- t(matrix(model$B/model$D, model$numGenes, length(predt)))

    Kxx <- multiKernCompute(simMultiKern, predt, t)
    predX <- c(meanPredX) + Re( Kxx%*%model$invK%*%model$m)
    varX <- Re( kernDiagCompute(simMultiKern, predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) ) 

    predExprs <- matrix(predX, length(predt), model$numGenes)
    varExprs <- matrix(varX, length(predt), model$numGenes)
    predExprs <- predExprs[1:(length(predt)-length(t)),]
    varExprs <- varExprs[1:(length(predt)-length(t)),]                
  }

  predF <- predF[1:(length(predt)-length(t))]
  varF <- varF[1:(length(predt)-length(t))]   

  model$predt <- predt[1:(length(predt)-length(t))]
  model$predF <- predF
  model$varF <- abs(varF)

  model$ypred <- predExprs
  model$ypredVar <- abs(varExprs)

  return (model)
}
