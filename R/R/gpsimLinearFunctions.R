gpsimCreate <- function(Ngenes, Ntf, times, y, yvar, options, genes = NULL) {

  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")
  
  if ( Ngenes != dim(y)[2] )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")
  
  model <- list(type="gpsim", y=as.array(y), yvar=as.array(yvar))

  kernType1 <- list(type="multi", comp=array())
  if ( any(grep("proteinPrior",names(options))) ) {
    model$proteinPrior <- options$proteinPrior
    kernType1$comp[1] <- "rbf"

    if ( any(grep("times",names(model$proteinPrior))) )
      timesCell <- list(protein=model$proteinPrior$times) else
    timesCell <- list(protein=times)

    tieParam <- 1

    for ( i in 1:Ngenes ) {
      kernType1$comp[i+1] <- "sim"
      timesCell <- c(timesCell, list(mRNA=times))
      tieParam <- c(tieParam, tieParam[length(tieParam)]+3)
    }

    model$timesCell <- timesCell
  } else {
    timesCell <- times
    tieParam <- 2

    kernType1 <- list(type="multi", comp=array())
    for ( i in 1:Ngenes ) {
      kernType1$comp[i] <- "sim"
      if ( i>1 ) tieParam <- c(tieParam, tieParam[length(tieParam)]+3)
    }
  }

  model$includeNoise <- options$includeNoise

  ## can be removed
  ## if ( model$includeNoise )
  ##   model$yvar <- array(0, length(yvar))

  if ( model$includeNoise ) {

    kernType2 <- list(type="multi", comp=array())
    if ( any(grep("proteinPrior",names(options))) ) {
      for ( i in 1:(Ngenes+1) )
        kernType2$comp[i] <- "white"
    } else {
      for ( i in 1:Ngenes )
        kernType2$comp[i] <- "white"        
    }     

    model$kern <- kernCreate(timesCell,
                             list(type="cmpnd", comp=list(kernType1, kernType2)))
  } else {
    model$kern <- kernCreate(timesCell, kernType1)
  }

  ## This is if we need to place priors on parameters ...
  ## ...
  ## ...

  model$kern <- modelTieParam(model$kern, tieParam)

  if ( any(grep("proteinPrior",names(options))) ) {
    for ( i in seq(2,model$kern$numBlocks,length.out=model$kern$numBlocks-1) ) {
      if ( model$includeNoise ) {
        model$D[i-1] <- model$kern$comp[[1]]$comp[[i]]$decay
        model$S[i-1] <- sqrt(model$kern$comp[[1]]$comp[[i]]$variance)
      } else {
        model$D[i-1] <- model$kern$comp[[i]]$decay
        model$S[i-1] <- sqrt(model$kern$comp[[i]]$variance)
      }
    }
  } else {
    for ( i in seq(length.out=model$kern$numBlocks) ) {
      if ( model$includeNoise ) {
        model$D[i] <- model$kern$comp[[1]]$comp[[i]]$decay
        model$S[i] <- sqrt(model$kern$comp[[1]]$comp[[i]]$variance)
      } else {
        model$D[i] <- model$kern$comp[[i]]$decay
        model$S[i] <- sqrt(model$kern$comp[[i]]$variance)
      }
    }
  }

  model$numParams <- Ngenes + model$kern$nParams
  model$numGenes <- Ngenes
  model$mu <- apply(y, 2, mean)
  model$B <- model$D*model$mu

  if ( any(grep("proteinPrior",names(model))) ) {
    model$m <- c(model$proteinPrior$values, model$y)
  } else {
    model$m <- c(model$y)
    model$t <- times
  }

  if ( any(grep("bTransform",names(options))) ) {
    model$bTransform <- options$bTransform
  } else {
    model$bTransform <- "positive"
  }

  model$optimiser <- options$optimiser

  if ( any(grep("fix", names(options))) ) {
    model$fix <- options$fix
  }

  if (!is.null(genes)) {
    model$genes <- genes
  }

  params <- gpsimExtractParam(model)
  model <- gpsimExpandParam(model, params)
  
  return (model)
 
}



gpsimExtractParam <- function (model, option=1) {
  ## option=1: only return parameter values;
  ## option=2: return both parameter values and names.

  funcName <- paste(optimiDefaultConstraint(model$bTransform), "Transform", sep="")
  func <- get(funcName, mode="function")

  if ( option == 1 ) {
    params <- kernExtractParam(model$kern)
    params <- c(params, func(model$B, "xtoa"))

    if ( any(grep("fix", names(model))) ) 
      for ( i in seq(along=model$fix$index) )
        params[model$fix$index[i]] <- model$fix$value[i]

    params <- Re(params)        

  } else {
    params <- kernExtractParam(model$kern, option)
    params$values <- c(params$values, func(model$B, "xtoa"))
    for ( i in seq(along=model$mu) ) {
      params$names <- c(params$names, paste("Basal", i, sep=""))
    }

    if ( any(grep("fix", names(model))) ) 
      for ( i in seq(along=model$fix$index) )
        params$values[model$fix$index[i]] <- model$fix$value[i]

    params$values <- Re(params$values)    
    
  }

  return (params)
}



gpsimExpandParam <- function (model, params) {
  if ( is.list(params) )
    params <- params$values

  params <- Re(params)
  if ( any(grep("fix", names(model))) ) 
    for ( i in seq(along=model$fix$index) )
      params[model$fix$index[i]] <- model$fix$value[i]

  if ( length(params) != model$numParams )
    stop("Parameter vector is incorrect length.")

  startVal <- 1
  endVal <- model$kern$nParams
  model$kern <- kernExpandParam(model$kern, params[startVal:endVal])

  funcName <- paste(optimiDefaultConstraint(model$bTransform), "Transform", sep="")
  func <- get(funcName, mode="function")

  model$B <- func(params[(endVal+1):length(params)], "atox")

  if ( any(grep("proteinPrior",names(model))) ) {
    for ( i in seq(2,model$kern$numBlocks,length.out=model$kern$numBlocks-1) ) {
      if ( model$includeNoise ) {
        model$D[i-1] <- model$kern$comp[[1]]$comp[[i]]$decay
        model$S[i-1] <- sqrt(model$kern$comp[[1]]$comp[[i]]$variance)
      } else {
        model$D[i-1] <- model$kern$comp[[i]]$decay
        model$S[i-1] <- sqrt(model$kern$comp[[i]]$variance)
      }
    }
  } else {
    for ( i in seq(length.out=model$kern$numBlocks) ) {
      if ( model$includeNoise ) {
        model$D[i] <- model$kern$comp[[1]]$comp[[i]]$decay
        model$S[i] <- sqrt(model$kern$comp[[1]]$comp[[i]]$variance)
      } else {
        model$D[i] <- model$kern$comp[[i]]$decay
        model$S[i] <- sqrt(model$kern$comp[[i]]$variance)
      }
    }
  }

  model$mu <- model$B/model$D

  model <- gpsimUpdateKernels(model)

  if ( any(grep("proteinPrior",names(model))) ) {
    if ( model$includeNoise ) {
      yInd <- 1:model$kern$comp[[1]]$diagBlockDim[2]
      mInd <- model$kern$comp[[1]]$diagBlockDim[1] + yInd

      for ( i in seq(length=model$numGenes) ) {
        model$m[mInd] <- model$y[yInd]-model$mu[i]*array(1, length(yInd), 1)
        yInd <- yInd+model$kern$comp[[1]]$diagBlockDim[i+1]
        mInd <- mInd+model$kern$comp[[1]]$diagBlockDim[i+1]
      }
    } else {
      yInd <- 1:model$kern$diagBlockDim[2]
      mInd <- model$kern$diagBlockDim[1] + yInd
      
      for ( i in seq(length=model$numGenes) ) {
        model$m[mInd] <- model$y[yInd]-model$mu[i]*array(1, length(yInd), 1)
        yInd <- yInd+model$kern$diagBlockDim[i+1]
        mInd <- mInd+model$kern$diagBlockDim[i+1]
      }         
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



gpsimUpdateKernels <- function (model) {
  eps <-  1e-6

  if ( any(grep("proteinPrior",names(model))) & any(grep("timesCell",names(model))) ) {
    k <- Re(kernCompute(model$kern, model$timesCell))
    if ( model$includeNoise ) {
      noiseVar <- c(array(0, model$kern$comp[[1]]$diagBlockDim[1], 1), model$yvar)
    } else {
      noiseVar <- c(array(eps, model$kern$comp[[1]]$diagBlockDim[1], 1), model$yvar)
    }
  } else {
    k <- Re(kernCompute(model$kern, model$t))
    noiseVar <- c(as.array(model$yvar))
  }

  model$K <- k+diag(as.array(noiseVar))
  invK <- jitCholInv(model$K)

  if ( is.nan(invK[1]) ) { 
    cat("kern$decay = \n", model$D, "\n")
    cat("\n")
    cat("kern$sensitivity = \n", model$S, "\n")
    cat("\n")
    cat("kern$flength = \n", model$kern$comp[1]$flength, "\n")
    cat("\n")
				
    stop("Singular chol(K) matrix!")
  }
  
  model$invK <- invK$invM

  if ( invK$jitter > 1e-4 )
    warning(paste("Warning: gpsimUpdateKernels added jitter of", signif(invK$jitter, digits=4)))

  model$logDetK <- 2* sum( log ( diag(invK$chol) ) )
 
  return (model)
}



cgpsimOptimise <- function (model, options, ...) {
  if ( any(grep("optimiser", names(options))) ) {
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
  if ( any(grep("bprior",names(model))) ) {
    ll <- ll + kernPriorLogProb(model$kern)
    ll <- ll + priorLogProb(model$bprior, model$B)
  }
  return (ll)
}



gpsimLogLikeGradients <- function (model) {
  covGrad <- -model$invK + model$invK %*% model$m %*% t(model$m) %*% model$invK
  covGrad <- 0.5*covGrad

  if ( any(grep("proteinPrior",names(model))) ) {
    g <- kernGradient(model$kern, model$timesCell, covGrad)
  } else {
    g <- kernGradient(model$kern, model$t, covGrad)
  }

  if ( any(grep("bprior",names(model))) ) {
    g <- g + kernPriorGradient(model$kern)
  }
  
  gmuFull <- t(model$m) %*% model$invK

  if ( any(grep("proteinPrior",names(model))) ) {
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

  funcName <- paste(optimiDefaultConstraint(model$bTransform), "Transform", sep="")
  func <- get(funcName, mode="function")

  ## prior contribution
  if ( any(grep("bprior",names(model))) ) {
    gb <- gb + priorGradient(model$bprior, model$B)
  }

  gb <- gb*func(model$B, "gradfact")

  ## gradient for D
  gd <- -gmu*model$B/(model$D*model$D)
  if ( model$kern$numBlocks == 1 ) {
    decayIndices <- 1
  } else if ( model$kern$numBlocks>1 ) {
    if ( any(grep("proteinPrior",names(model))) ) {
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

  if ( any(grep("fix",names(model))) ) 
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



cgpsimExtractParam <- function (model, option=1) {
  ## option=1: only return parameter values;
  ## option=2: return both parameter values and names.

  if ( option == 1 ) {
    params <- gpsimExtractParam(model$comp[[1]])
  } else {
    params <- gpsimExtractParam(model$comp[[1]], option)
  }

  return (params)
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



gpsimGradient <- function (params, model, ...) {
  model <- gpsimExpandParam(model, params)
  g <- - gpsimLogLikeGradients (model)

  return (g)
}



gpsimUpdateProcesses <- function (model) {
  if ( any(grep("proteinPrior",names(model))) ) {
    t <- model$timesCell[[2]]
  } else {
    t <- model$t
  }

  numData <- length(t)
  predt <- c(seq(min(t), max(t), length=100), t)
  ymean <- t(matrix(model$y[1,], model$numGenes, numData))

  if ( any(grep("proteinPrior",names(model))) ) {
    predTimeCell <- list()
    if ( model$includeNoise ) {
      for ( i in seq(length=model$kern$comp[[1]]$numBlocks) )
        predTimeCell[[i]] <- predt

      Kxx <- multiKernCompute(model$kern$comp[[1]], predTimeCell, model$timesCell)
      diagKxx <- kernDiagCompute(model$kern$comp[[1]], predTimeCell)
      x <- c(model$proteinPrior$values, model$y-ymean)
      ind <- 1:length(predTimeCell[[1]])
      predFull <- list()
      varFull <- list()
      for ( indBlock in seq(length=model$kern$comp[[1]]$numBlocks) ) {
        K <- Kxx[ind,]
        diagK <- diagKxx[ind,]
        predFull[[indBlock]] <- Re( K %*% model$invK %*% x )
        varFull[[indBlock]] <- Re( diagK - t(apply(t(K)*(model$invK %*% t(K)), 2, sum)) )
        ind <- ind + length(predTimeCell[[indBlock]])
      }
    } else {
      for ( i in seq(length=model$kern$numBlocks) )
        predTimeCell[[i]] <- predt

      Kxx <- multiKernCompute(model$kern, predTimeCell, model$timesCell)
      diagKxx <- kernDiagCompute(model$kern, predTimeCell)
      x <- c(model$proteinPrior$values, model$y-ymean)

      ind <- 1:length(predTimeCell[[1]])
      predFull <- list()
      varFull <- list()
      for ( indBlock in seq(length=model$kern$numBlocks) ) {
        K <- Kxx[ind,]
        diagK <- diagKxx[ind,]
        predFull[[indBlock]] <- Re( K %*% model$invK %*% x )
        varFull[[indBlock]] <- Re( diagK - t(apply(t(K)*(model$invK%*%t(K)), 2)) )
        ind <- ind + length(predTimeCell[[indBlock]])
      }
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
    if ( model$includeNoise ) {
      proteinKern$inverseWidth <- model$kern$comp[[1]]$comp[[1]]$inverseWidth
    } else {
      proteinKern$inverseWidth <- model$kern$comp[[1]]$inverseWidth
    }

    if ( model$includeNoise) {
      K <- simXrbfKernCompute(model$kern$comp[[1]]$comp[[1]], proteinKern, model$t, predt)
      for ( i in seq(2, length=(model$kern$numBlocks-1)) ) 
        K <- rbind(K, simXrbfKernCompute(model$kern$comp[[1]]$comp[[i]], proteinKern, model$t, predt))
    } else {
      K <- simXrbfKernCompute(model$kern$comp[[1]], proteinKern, model$t, predt)
      for ( i in seq(2, length=(model$kern$numBlocks-1)) ) 
        K <- rbind(K,simXrbfKernCompute(model$kern$comp[[i]], proteinKern, model$t, predt))    }

    predF <- Re( t(K)%*%model$invK%*%c(model$y-ymean) )
    varF <- Re( kernDiagCompute(proteinKern, predt) - apply(K*(model$invK%*%K), 2, sum) )
    meanPredX <- t(matrix(model$B/model$D, model$numGenes, length(predt)))

    if ( model$includeNoise ) {
      Kxx <- multiKernCompute(model$kern$comp[[1]], predt, t)
      predX <- c(meanPredX) + Re( Kxx%*%model$invK%*%c(model$y-ymean) )
      varX <- Re( kernDiagCompute(model$kern$comp[[1]], predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) ) 
    } else {
      Kxx <- multiKernCompute(model$kern, predt, t)
      predX <- c(meanPredX) + Re( Kxx%*%model$invK%*%c(model$y-ymean) )
      varX <- Re( kernDiagCompute(model$kern, predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) ) 
    }

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



gpsimBarencoResults <- function (model, scale, expType, expNo, option=1) {
  geneNames <- c("DDB2", "hPA26", "TNFRSF10b", "p21", "BIK")
  order <- c(1,5,3,4,2)

  modelB <- model$comp[[1]]$B*scale
  scaleModelB <- modelB/mean(modelB)
  modelS <- model$comp[[1]]$S*scale
  modelS <- modelS/modelS[4]
  
  ## Display the result
  if ( option==1 ) {
    for ( i in seq(along=model$comp) ) {
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, ylim=c(-1,3), type="l", lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

      for ( j in seq(length=model$comp[[i]]$numGenes) ) {
        x11()
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(0,4), type="l", lwd=3, xlab="Time",ylab="")
##        title(paste("mRNA", geneNames[order[j]]))
        points(model$comp[[i]]$timesCell[[2]], model$comp[[i]]$y[,j], lwd=3, col=3) 
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
    }

    x11()
    barplot(modelB[order], names.arg=geneNames, mfg=c(3,1))
    title("Basal Transcription Rates")
    x11()
    barplot(model$comp[[1]]$D[order], names.arg=geneNames, mfg=c(2,1))
    title("Decay rate of the mRNA")
    x11()
    barplot(modelS[order], names.arg=geneNames, mfg=c(1,1))
    title("Sensitivities to the Transcription Factor")
    
  } else {
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
    for ( i in seq(along=model$comp) ) {
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, ylim=c(-1,3), type="l", lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53 using a Linear Response Model")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

      for ( j in seq(length=model$comp[[i]]$numGenes) ) {
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(0,4), type="l", lwd=3, xlab="Time",ylab="")
        title(paste("mRNA", geneNames[order[j]]))
        points(model$comp[[i]]$timesCell[[2]], model$comp[[i]]$y[,j], lwd=3, col=3)       
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }

    }

    barplot(modelB[order], names.arg=geneNames, mfg=c(3,1))
    title("Basal Transcription Rates")
    barplot(model$comp[[1]]$D[order], names.arg=geneNames, mfg=c(2,1))
    title("Decay rate of the mRNA")
    barplot(modelS[order], names.arg=geneNames, mfg=c(1,1))
    title("Sensitivities to the Transcription Factor")
    
    dev.off()
  }

}



gpsimElkResults <- function (model, scale, genes, expType, expNo, option=1) {
  ## Display the result
  if ( option==1 ) {
    for ( i in seq(along=model$comp) ) {
      ymin <- min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF))
      ymax <- max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, type="l", lwd=3, ylim=c(ymin, ymax), xlab="Time",ylab="")
      title("Predicted Protein Concentration")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      
      for ( j in seq(length=model$comp[[i]]$numGenes) ) {
        x11()
        ymin <- min(c(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), model$comp[[i]]$y[,j]))
        ymax <- max(c(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), model$comp[[i]]$y[,j]))
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], type="l", lwd=3, ylim=c(ymin, ymax), xlab="Time",ylab="")
        title(paste("mRNA", genes[j]))
        if ( any(grep("proteinPrior",names(options))) ) {
          points(model$comp[[i]]$timesCell[[2]], model$comp[[i]]$y[,j], lwd=3, col=3)
        } else {
          points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
        }
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
    }

    x11()
    barplot(model$comp[[1]]$B, mfg=c(3,1))
    title("Basal Transcription Rates")
    x11()
    barplot(model$comp[[1]]$D, mfg=c(2,1))
    title("Decay rate of the mRNA")
    x11()
    barplot(model$comp[[1]]$S, mfg=c(1,1))
    title("Sensitivities to the Transcription Factor")
    
  } else {
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
    for ( i in seq(along=model$comp) ) {
      ymin <- min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF))
      ymax <- max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))
      
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, type="l", lwd=3, ylim=c(ymin, ymax), xlab="Time",ylab="")
      title("Predicted Protein Concentration using a Linear Response Model")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      
      for ( j in seq(length=model$comp[[i]]$numGenes) ) {
        ymin <- min(c(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), model$comp[[i]]$y[,j]))
        ymax <- max(c(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), model$comp[[i]]$y[,j]))
     
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], type="l", lwd=3, ylim=c(ymin, ymax), xlab="Time",ylab="")
        title(paste("mRNA", genes[j]))
        if ( any(grep("proteinPrior",names(options))) ) {
          points(model$comp[[i]]$timesCell[[2]], model$comp[[i]]$y[,j], lwd=3, col=3)
        } else {
          points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
        }
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }

    }

    barplot(model$comp[[1]]$B, mfg=c(3,1))
    title("Basal Transcription Rates")

    barplot(model$comp[[1]]$D, mfg=c(2,1))
    title("Decay rate of the mRNA")

    barplot(model$comp[[1]]$S, mfg=c(1,1))
    title("Sensitivities to the Transcription Factor")
    
    dev.off()
  }

}
