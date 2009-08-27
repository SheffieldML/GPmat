gpdisimCreate <- function(Ngenes, Ntf, times, y, yvar, options, annotation=NULL) {

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



gpdisimExtractParam <- function (model, option=1) {
  params <- gpsimExtractParam(model, option)
  return (params)
}



gpdisimExpandParam <- function (model, params) {

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

  model$delta <- model$kern$comp[[2]]$di_decay
  model$sigma <- model$kern$comp[[2]]$di_variance

  for ( i in seq(2, model$kern$numBlocks) ) {
    model$D[i-1] <- model$kern$comp[[i]]$decay
    model$S[i-1] <- sqrt(model$kern$comp[[i]]$variance)
  }

  model$mu <- model$B/model$D
  model <- gpsimUpdateKernels(model)

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
  if ( any(grep("bprior",names(model))) ) {
    ll <- ll + kernPriorLogProb(model$kern)
    ll <- ll + priorLogProb(model$bprior, model$B)
  }
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
    ind <- ind + numData
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
  decayIndices <- 5

  for ( i in seq(3, model$kern$numBlocks, length=(model$kern$numBlocks-2)) )
    decayIndices <- c(decayIndices, decayIndices[length(decayIndices)]+2)

  g[decayIndices] <- g[decayIndices]+gd*expTransform(model$D, "gradfact")

  g <- c(g, gb)

  if ( any(grep("fix",names(model))) )
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



cgpdisimExtractParam <- function (model, option=1) {
  ## option=1: only return parameter values;
  ## option=2: return both parameter values and names.

  if ( option == 1 ) {
    params <- gpdisimExtractParam(model$comp[[1]])
  } else {
    params <- gpdisimExtractParam(model$comp[[1]], option)
  }

  return (params)
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



gpdisimUpdateProcesses <- function (model) {

  if ( any(grep("proteinPrior",names(model))) ) {
    t <- model$timesCell[[2]]
  } else {
    t <- model$t
  }

  numData <- length(t)
  predt <- c(seq(min(t), max(t), length=100), t)
  proteinKern <- kernCreate(model$t, "sim")
  proteinKern$inverseWidth <- model$kern$comp[[1]]$inverseWidth
  proteinKern$decay <- model$delta
  proteinKern$variance <- model$kern$comp[[2]]$di_variance
  K <- simXrbfKernCompute(proteinKern, model$kern$comp[[1]], predt, model$t)
  for ( i in seq(2, length=(model$kern$numBlocks-1)) )
    K <- cbind(K,t(disimXsimKernCompute(model$kern$comp[[i]], proteinKern, model$t, predt)))

  ymean <- t(matrix(c(0, model$mu), model$numGenes+1, numData))

  predF <- K %*% model$invK %*% c(model$y-ymean)
  #varF <- kernDiagCompute(proteinKern, predt) - apply(t(K)*(model$invK %*% t(K)), 2, sum)
  varF <- model$kern$comp[[1]]$variance * kernDiagCompute(proteinKern, predt) - apply(t(K)*(model$invK %*% t(K)), 2, sum)

  Kxx <- multiKernCompute(model$kern, predt, model$t)
  meanPredX <- t(matrix(c(0, model$B/model$D), model$numGenes+1, length(predt)))
  predX <- c(meanPredX) + Re(Kxx %*% model$invK %*% c(model$y-ymean))
  varX <- Re( kernDiagCompute(model$kern, predt)-apply(t(Kxx)*(model$invK%*%t(Kxx)), 2, sum) )

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

  model$predt <- predt
  model$predF <- predF
  model$varF <- abs(varF)

  model$ypred <- predExprs
  model$ypredVar <- abs(varExprs)

  return (model)
}






gpdisimMef2Results <- function (model, scale, expType, expNo, option=1) {
  geneNames <- c("gene1", "Rya-r44F", "gene3 ", "ttk", "gene5", "gene6")

  ## Display the result
  if ( option==1 ) {
    for ( i in seq(along=model$comp) ) {
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, type="l", lwd=3, xlab="Time",ylab="")
      title("Inferred Mef2 Protein Concentration")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

      for ( j in seq(length=model$comp[[i]]$numGenes+1) ) {
        x11()
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(0,4), type="l", lwd=3, xlab="Time",ylab="")
        if ( j==1 ) {
          title(paste("Driving Input mRNA"))
        } else {
          title(paste("mRNA", geneNames[j]))
        }
        points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
    }

  } else {
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
    for ( i in seq(along=model$comp) ) {
      plot(model$comp[[i]]$predt, model$comp[[i]]$predF, type="l", lwd=3, xlab="Time",ylab="")
      title("Inferred Mef2 Protein Concentration")
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

      for ( j in seq(length=model$comp[[i]]$numGenes) ) {
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(0,4), type="l", lwd=3, xlab="Time",ylab="")
        if ( j==1 ) {
          title(paste("Driving Input mRNA"))
        } else {
          title(paste("mRNA", geneNames[j]))
        }
        points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
    }

    dev.off()
  }

}

