gpsimMapCreate <- function (numGenes, numProteins, times, y, yvar, options) {
  if ( any(dim(y)!=dim(yvar)) )
    stop("The gene variances have a different size matrix to the gene values.")
  
  if ( numGenes != dim(y)[2] )
    stop("The number of genes given does not match the dimension of the gene values given.")

  if ( length(times) != dim(y)[1] )
    stop("The number of time points given does not match the number of gene values given.")

  if ( length(options$B) != numGenes )
    stop("Incorrect length of options$B.")
  if ( length(options$D) != numGenes )
    stop("Incorrect length of options$D.")
  if ( length(options$S) != numGenes )
    stop("Incorrect length of options$S.")

  model <- list(type="gpsimMap", y=as.array(y), yvar=as.array(yvar), t=times, nonLinearity=options$nonLinearity, Transform="exp", B=options$B, D=options$D, S=options$S, optimiser=options$optimiser)

  if ( any(grep("gParam", names(options))) ) {
    model$gParam=options$gParam
    model$ngParam <- dim(as.matrix(options$gParam))[1]*numGenes
  } else {
    model$ngParam <- 0
  }

  if ( any(grep("bTransform", names(options))) )
    model$bTransform <- options$bTransform
  if ( any(grep("alphaTransform", names(options))) )
    model$alphaTransform <- options$alphaTransform

  if ( any(names(model)=="gParam") ) {
    model$nonLinearity <- array()
    if ( !is.list(options$nonLinearity) ) {
      for ( i in seq(length=numGenes) ) {
        model$nonLinearity[i] <- options$nonLinearity
        model$isGroupNonlinearity <- 1
      }
    } else {
      model$isGroupNonlinearity <- 0
      model$nonLinearity <- options$nonLinearity
    }
  }

  epsilon <- 1e-2

  if ( any(names(options)=="alpha") ) {
    model$alpha <- options$alpha
    model$includeRepression <- 1
  } else if ( model$nonLinearity!="linear" ) {
    if ( any(grep("isGroupNonlinearity", names(model))) ) {
      if ( model$isGroupNonlinearity )
        if ( model$nonLinearity[1]=="repression" ) {
          model$alpha <- epsilon*array(1, numGenes)
          model$includeRepression <- 1
        } else {
          model$alpha <- array(0, numGenes)
          model$includeRepression <- 0      
        }
    } else {
      for ( i in seq(length=numGenes) )
        if ( model$nonLinearity[i]=="repression" ) {
          model$alpha[i] <- epsilon
          model$includeRepression <- 1       
        } else {
          model$alpha[i] <- 0
        }
    }
  }

  noiseInit <- 1e-2
  if ( options$includeNoise ) {
 #   model$yvar <- array(0, dim=dim(model$yvar))
    model$includeNoise <- options$includeNoise
    model$noiseVar <- noiseInit*array(1, numGenes)
  }

  if ( any(grep("startPoints", names(options))) ) {
    startPoint <- options$startPoint
    endPoint <- options$endPoint
    span <- endPoint - startPoint
    fullSpan <- endPoint - startPoint
    if ( any(grep("intPoints", names(options))) ) {
      model$step <- fullSpan/(options$intPoints-1)
    } else { 
      model$step <- options$step
    }

    model$mapt <- seq(startPoint, times[1], by=model$step)

    model$timesIndex <- length(model$mapt)
    for ( i in seq(length=(length(times)-1)) ) {
      model$mapt <- c(model$mapt, seq(times[i]+model$step, times[i+1], by=model$step))
      model$timesIndex <- c(model$timesIndex, length(model$mapt))
    }
    
    if ( endPoint > times[length(times)] )
      model$mapt <- c(model$mapt, seq(times[length(times)]+model$step, endPoint, by=model$step))
    
  } else if ( any(grep("mapt", names(options))) ) {
    model$mapt <- options$mapt
    model$step <- options$step
    model$timesIndex <- options$timesIndex
  }

  model$numMapPts <- length(model$mapt)
  model$numGenes <- numGenes  

  if ( any(grep("fix", names(options))) ) {
    model$fix <- options$fix
  }
  
  if ( any(grep("proteinPrior",names(options))) ) {
    model$proteinPrior <- options$proteinPrior
    model$consLambda <- 0    
  }

  if ( is.list(options$kern) ) {
    model$kern <- options$kern
  } else {
    if ( length(options$kern) > 1 ) {
      kernType <- list(type=options$kern[1])
      for ( i in seq(length=(length(options$kern)-1)) )
        kernType$comp[i] <- options$kern[i+1]
    } else {
      kernType <- options$kern
    }
    model$kern <- kernCreate(model$mapt, kernType)
  }

  model$f <- array(0, length(model$mapt))
  model$varf <- array(1, length(model$mapt))

  model <- gpsimMapUpdateG(model)
  model <- gpsimMapUpdateYpred(model)

  model$W <- matrix(0, length(model$mapt), length(model$mapt))
  model$updateW <- TRUE

  model <- gpsimMapUpdateKernels(model)
  model <- gpsimMapFunctionalUpdateW(model)

  param <- gpsimMapExtractParam(model)
  model$numParams <- length(param)

  return (model)
}



gpsimMapUpdateG <- function (model) {
  f <- as.array(model$f)

  if ( model$ngParam == 0 ) {
    if ( model$nonLinearity == "linear" ) {
      model$g <- f
      model$gGrad <- array(1, length(model$g))
      model$gGrad2 <- array(0, length(model$g))
      model$gGrad3 <- array(0, length(model$g))
      model$isConcave <- TRUE
    } else if ( model$nonLinearity == "exp" ) {
      model$g <- exp(f)
      model$gGrad <- model$g
      model$gGrad2 <- model$g
      model$gGrad3 <- model$g
      model$isConcave <- TRUE
    } else if ( model$nonLinearity == "quadratic" ) {
      model$g <- f*f
      model$gGrad <- 2*f
      model$gGrad2 <- 2*array(1, length(model$g))
      model$isConcave <- FALSE
    } else if ( model$nonLinearity == "negLogLogit" ) {
      model$g <- log(1+exp(f))
      model$gGrad <- sigmoid(f)
      model$gGrad2 <- model$gGrad*(1-model$gGrad)
      model$isConcave <- TRUE
    } else if ( model$nonLinearity == "sigmoid" ) {
      model$g <- sigmoid(f)
      model$gGrad <- model$g*(1-model$g)
      model$gGrad2 <- model$g*(1-model$g) - 2*model$g*(model$g*(1-model$g))
      model$isConcave <- FALSE
    } else if ( model$nonLinearity == "repression" ) {
      model$g <- gpsimModelFunctions(model,"g")
      model$gGrad <- gpsimModelFunctions(model,"grad")
      model$gGrad2 <- gpsimModelFunctions(model,"grad2")
      model$gGrad3 <- gpsimModelFunctions(model,"grad3")   
      model$isConcave <- TRUE
    } else if ( model$nonLinearity == "activation" ) {
      model$g <- gpsimModelFunctions(model,"g")
      model$gGrad <- gpsimModelFunctions(model,"grad")
      model$gGrad2 <- gpsimModelFunctions(model,"grad2")
      model$gGrad3 <- gpsimModelFunctions(model,"grad3")   
      model$isConcave <- TRUE
      if ( !any(is.na(model$gParam)) ) {
        model$ngParam <- length(model$gParam[1,])
        model$dg <- gpsimModelFunctions(model,"paramGrad")
        model$dg2 <- gpsimModelFunctions(model,"paramGrad2")
      }
    } else {
      stop("Invalid non-linearity.")
    }
  } else {
    model$g <- gpsimModelFunctions(model,"g")
    model$gGrad <- gpsimModelFunctions(model,"grad")
    model$gGrad2 <- gpsimModelFunctions(model,"grad2")
    model$gGrad3 <- gpsimModelFunctions(model,"grad3")   
    model$isConcave <- TRUE
    model$dg <- gpsimModelFunctions(model,"paramGrad")
    model$dggrad <- gpsimModelFunctions(model,"paramGgrad")
    model$dggrad2 <- gpsimModelFunctions(model,"paramGgrad2")
  
  }
 
  return (model)
}



gpsimModelFunctions <- function (model, options="g") {

  if ( model$ngParam ) {
    gf <- matrix(0, length(model$f), model$numGenes)
    grad <- matrix(0, length(model$f), model$numGenes)
    grad2 <- matrix(0, length(model$f), model$numGenes)

    for ( k in seq(length=model$numGenes) ) {
      gParamk <- as.matrix(model$gParam)[,k]
      if ( model$nonLinearity[k]=="repression" ) {
        gamma <- gParamk[1]
        
        expf <- expTransform(model$f, "atox")
        exp2f <- expf*expf
        exp3f <- expf*exp2f        

        if ( options=="g" ) {
          gf[,k] <- 1/(gamma+expf)
        } else if ( options=="grad" ) {
          gf[,k] <- -expf/((gamma+expf)*(gamma+expf))
        } else if ( options=="grad2" ) {
          grad[,k] <- -expf/((gamma+expf)*(gamma+expf))
          gf[,k] <- grad[,k] + 2*exp2f/(gamma+expf)^3
        } else if ( options=="grad3" ) {
          grad[,k] <- -expf/((gamma+expf)*(gamma+expf))
          grad2[,k] <- grad[,k] + 2*exp2f/(gamma+expf)^3
          gf[,k] <- grad2[,k] + 4*exp2f/(gamma+expf)^3 - 6*exp3f/(gamma+expf)^4
        } else if ( options=="paramGrad" ) {
          gf[,k] <- -1/(gamma+expf)^2
        } else if ( options=="paramGgrad" ) {
          gf[,k] <- 2*expf/((gamma+expf)^3)
        } else if ( options=="paramGgrad2" ) {
          gf[,k] <- 2*expf/((gamma+expf)^3) - 6*exp2f/((gamma+expf)^4)
        }
      } else if ( model$nonLinearity[k]=="activation" ) {     
        gamma <- gParamk[1]
         
        expmf <- expTransform(-model$f, "atox")
        expm2f <- expmf*expmf
        expm3f <- expmf*expm2f

        if ( options=="g" ) {
          gf[,k] <- 1/(1+gamma*expmf)
        } else if ( options=="grad" ) {
          gf[,k] <- gamma*expmf/((1+gamma*expmf)*(1+gamma*expmf))
        } else if ( options=="grad2" ) {
          grad[,k] <- -gamma*expmf/((1+gamma*expmf)*(1+gamma*expmf))
          gf[,k] <- grad[,k] + 2*gamma^2*expm2f/(1+gamma*expmf)^3
        } else if ( options=="grad3" ) {
          grad[,k] <- -gamma*expmf/((1+gamma*expmf)*(1+gamma*expmf))
          grad2[,k] <- grad[,k] + 2*gamma^2*expm2f/(1+gamma*expmf)^3
          gf[,k] <- -grad2[,k] - 4*gamma^2*expm2f/(1+gamma*expmf)^3 + 6*gamma^3*expm3f/(1+gamma*expmf)^4
        } else if ( options=="paramGrad" ) {
          gf[,k] <- -expmf/((1+gamma*expmf)*(1+gamma*expmf))
        } else if ( options=="paramGgrad" ) {
          grad <- expmf/((1+gamma*expmf)^2)
          gf[,k] <- grad - 2*gamma*expm2f/((1+gamma*expmf)^3)
        } else if ( options=="paramGgrad2" ) {
          gf[,k] <- expmf/((1+gamma*expmf)^2)
          gf[,k] <- gf[,k] - 2*gamma*expm2f/((1+gamma*expmf)^3)
          gf[,k] <- -gf[,k] + 4*gamma*expm2f/((1+gamma*expmf)^3) - 6*gamma^2*expm3f/((1+gamma*expmf)^4)            
        }
      }
    }
  }

  return (gf)
}



gpsimMapUpdateYpred <- function (model, option=1) {
  ## option=1: normal integral computation
  ## option=2: Simpson's rule integral computation

  start <- 1
  model$ypred <- matrix(0, model$numMapPts, model$numGenes)
  for ( j in seq(length=model$numGenes) ) {
    lnintegral <- array(0, model$numMapPts)
    weight <- array(0, model$numMapPts)
    
    if ( any(grep("ngParam", names(model))) ) {
      if ( model$ngParam ) {
        gInd <- j
      } else {
        gInd <- 1
      }
    } else {
      gInd <- 1
    }

    for ( i in seq(length=model$numMapPts) ) {
      tfs <- model$mapt[i]      
      if ( any(grep("isGroupNonlinearity", names(model))) ) {
        base <- model$alpha[j]*exp(-model$D[j]*tfs) 
      } else {
        base <- 0
      }

      if ( option==1 ) {
        lnintegral[i] <- (model$D[j]*tfs)+log(model$step)
        model$ypred[i,j] <- base + model$B[j]/model$D[j]+sum(as.matrix(model$g)[start:i, gInd]*exp(log(model$S[j])-model$D[j]*tfs+lnintegral[start:i]))        

      } else if ( option==2 ) {
        
        lnintegral[i] <- model$D[j]*tfs + log(model$step)
        weight[i] <- 1
        if ( i%%2 ) {          
          w <- weight
          w[seq(2, length(w)-1, by=2)] <- 4
          w[seq(3, length(w)-2, by=2)] <- 2
        } else {
          w <- weight
          w[1] <- 3
          w[seq(3, length(w)-1, by=2)] <- 4
          w[seq(4, length(w)-2, by=2)] <- 2
        }
        
        model$ypred[i,j] <- base + model$B[j]/model$D[j] + sum(model$g[start:i, gInd]*w[start:i]*exp(log(model$S[j])-model$D[j]*tfs+lnintegral[start:i]))/3
      }
    }
  }
  
  return (model)
}



gpsimMapUpdateKernels <- function (model) {

  model$K <- kernCompute(model$kern, model$mapt)
  model$K <- model$K + 1e-6*diag(dim(model$K)[1])

  invK <- jitCholInv(model$K)
  model$invK <- invK$invM

  if ( invK$jitter > 1e-4 )
    warning(paste("Warning: gpsimUpdateKernels added jitter of", signif(invK$jitter, digits=4)))

  model$logDetK <- 2* sum( log ( diag(invK$chol) ) )

  model <- gpsimMapUpdatePosteriorCovariance(model)
  
  return (model)
}



gpsimMapUpdatePosteriorCovariance <- function (model) {
  ev <- eigen(model$W)

  U <- ev$vectors
  Lambda <- ev$values

  sortLambda <- sort(Lambda, decreasing=TRUE, index.return=TRUE)
  lambda <- sortLambda$x
  order <- sortLambda$ix

  U <- U[,order]

  ## deal with non positive definite matrix
  lambda[grep(TRUE, lambda<0)] <- 0
  ind <- seq(length=length(grep(TRUE, lambda!=0)))
  LambdaHalf <- diag(sqrt(lambda))
  LambdaHalf <- LambdaHalf[ind, ind]
  if ( length(ind)==1 ) {
    model$Whalf <- as.matrix(U[,ind]*LambdaHalf)
  } else {
    model$Whalf <- U[,ind]%*%LambdaHalf
  }

  if ( dim(model$Whalf)[2] == 0 )
    model$Whalf <- matrix(0, dim(model$Whalf)[1], 1)

  model$invCovf <- model$invK + model$Whalf %*% t(model$Whalf)
  KWhalf <- model$K %*% model$Whalf
  inner <- diag(dim(model$Whalf)[2]) + t(model$Whalf) %*% KWhalf

  innerInvSet <- jitCholInv(inner)
  innerInv <- innerInvSet$invM
  logDetInner <- 2* sum( log ( diag(innerInvSet$chol) ) )

  model$covf <- model$K - KWhalf %*% innerInv %*% t(KWhalf)
  model$logDetCovf <- model$logDetK - logDetInner

  if ( any(grep("proteinPrior",names(options))) ) {
    nCons <- length(model$proteinPrior$values)
    consMat <- as.matrix(array(0, dim(model$covf)))
    for ( k in seq(length=nCons) ) {
      ftimeIndex <- grep(TRUE, (model$proteinPrior$times[k]-model$mapt)==0)
      consMat[ftimeIndex, ftimeIndex] <- model$consLambda
    }
    hf <- solve(model$invCovf+consMat)
  } else {
    hf <- model$covf
  }

  model$varf <- diag(hf)
  
  return (model)
}



gpsimMapFunctionalUpdateW <- function (model) {
  numData <- length(model$t)
  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
  } else {
    yvar <- model$yvar
  }

  if ( model$updateW ) {
    intPoints <- (model$timesIndex[1]+1):model$numMapPts

    step2 <- model$step*model$step
    S2 <- model$S*model$S
    for ( k in intPoints ) {
      for ( l in intPoints ) {
        temp <- 0
        for ( i in seq(length=numData) ) {
          arg1 <- model$t[i]-model$mapt[k]
          arg2 <- model$t[i]-model$mapt[l]
          if ( arg1>=0 & arg2>=0 )
            for ( j in seq(length=model$numGenes) ) {
              if ( model$ngParam ) {
                gInd <- j
              } else {
                gInd <- 1
              }
              ind <- i + (j-1)*numData
              betaij <- 1/yvar[ind]
              temp <- temp+betaij*as.matrix(model$gGrad)[k,gInd]*as.matrix(model$gGrad)[l,gInd]*exp(-model$D[j]*(arg1+arg2)+log(S2[j])+log(step2))
            }
        }
        model$W[k,l] <- temp

        if ( is.na(model$W[k,l]) )
          warning("model$W is NA.")
      }
    }

    for ( k in intPoints ) {
      temp <- 0
      for ( i in seq(length=numData) ) {
        arg <- model$t[i]-model$mapt[k]
        if ( arg>=0 )
          for ( j in seq(length=model$numGenes) ) {
            if ( model$ngParam ) {
              gInd <- j
            } else {
              gInd <- 1
            }
            ind <- i+(j-1)*numData
            betaij <- 1/yvar[ind]
            factor <- (model$ypred[model$timesIndex[i], j]-model$y[ind])*betaij
            temp <- temp + factor*as.matrix(model$gGrad2)[k,gInd]*exp(-model$D[j]*arg+log(model$S[j])+log(model$step))
          }
      }
      model$W[k,k] <- model$W[k,k]+temp
      
      if ( is.na(model$W[k,k]) )
        warning("model$W is NA.")
    }

    model$updateW <- TRUE
    model$W <- 0.5*(model$W+t(model$W))
  }

  if ( all(model$nonLinearity == "linear") )
    model$updateW <- FALSE
  
  return (model)
}



gpsimMapExtractParam <- function (model, only.values=TRUE) {
  funcName <- paste(model$Transform, "Transform", sep="")
  func <- get(funcName, mode="function")

  if ( only.values ) {
    params <- kernExtractParam(model$kern)
    for ( i in seq(length=model$numGenes) ) {
      if ( any(grep("bTransform", names(model))) ) {
        if ( is.na(model$bTransform) ) {
          params <- c(params, model$B[i])
        } else {
          params <- c(params, func(model$B[i], "xtoa"))
        }
      } else {
        params <- c(params, func(model$B[i], "xtoa"))
      }

      params <- c(params, func(model$S[i], "xtoa"), func(model$D[i], "xtoa"))

      if ( any(grep("includeRepression", names(model))) )
        if ( model$includeRepression ) {
          if ( any(grep("alphaTransform", names(model))) ) {
            if ( is.na(model$alphaTransform) ) {
              params <- c(params, model$alpha[i])
            } else {
              params <- c(params, func(model$alpha[i], "xtoa"))
            }
          } else {
            params <- c(params, func(model$alpha[i], "xtoa"))
          }
        }

      if ( model$ngParam > 0 )
        params <- c(params, func(as.matrix(model$gParam)[,i], "xtoa"))

    }

    if ( model$includeNoise )
      params <- c(params, func(sqrt(model$noiseVar), "xtoa"))

    if ( any(grep("fix", names(model))) ) 
      for ( i in seq(along=model$fix$index) )
        params[model$fix$index[i]] <- model$fix$value[i]

    params <- Re(params)
    
  } else {
    params <- kernExtractParam(model$kern, only.values=FALSE)
    for ( i in seq(length=model$numGenes) ) {
      if ( any(grep("bTransform", names(model))) ) {
        if ( is.na(model$bTransform) ) {
          params$values <- c(params$values, model$B[i])
        } else {
          params$values <- c(params$values, func(model$B[i], "xtoa"))
        }
      } else {
        params$values <- c(params$values, func(model$B[i], "xtoa"))
      }

      params$values <- c(params$values, func(model$S[i], "xtoa"), func(model$D[i], "xtoa"))

      if ( any(grep("includeRepression", names(model))) )
        if ( model$includeRepression ) {
          if ( any(grep("alphaTransform", names(model))) ) {
            if ( is.na(model$alphaTransform) ) {
              params$values <- c(params$values, model$alpha[i])
            } else {
              params$values <- c(params$values, func(model$alpha[i], "xtoa"))
            }
          } else {
            params$values <- c(params$values, func(model$alpha[i], "xtoa"))
          }
        }

      params$names <- c(params$names, paste("Basal", i, sep=""))
      params$names <- c(params$names, paste("Sensitivity", i, sep=""))
      params$names <- c(params$names, paste("Decay", i, sep=""))

      if ( any(grep("isGroupNonlinearity", names(model))) ) {
        if ( model$isGroupNonlinearity ) {
          if ( model$nonLinearity[i] == "repression" )
            params$names <- c(params$names, paste("alpha", i, sep=""))
        } else {
          params$names <- c(params$names, paste("alpha", i, sep=""))
        }
      } else {
        params$names <- c(params$names, paste("alpha", i, sep=""))
      }
      
      if ( model$ngParam > 0 ) {
        params$values <- c(params$values, func(as.matrix(model$gParam)[,i], "xtoa"))
        for ( k in seq(along=as.matrix(model$gParam)[,i]) )
          params$names <- c(params$names, paste("gParams", i, sep=""))
      }
    }

    if ( model$includeNoise ) {
      params$values <- c(params$values, func(sqrt(model$noiseVar), "xtoa"))
      for ( i in seq(length=model$numGenes) )
        params$names <- c(params$names, paste("noiseSd", i, sep=""))
    }

    if ( any(grep("fix", names(model))) ) 
      for ( i in seq(along=model$fix$index) )
        params$values[model$fix$index[i]] <- model$fix$value[i]

    params$values <- Re(params$values)
    trueparams <- params$values
    names(trueparans) <- params$names
    params <- trueparams
  }

  return (params)
}



gpsimMapExpandParam <- function (model, params) {

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

  ngParamk <- model$ngParam/model$numGenes

  if ( !any(grep("isGroupNonlinearity", names(model))) ) {
    nBaseParam <- 3
  } else {
    if ( model$isGroupNonlinearity & (model$nonLinearity[1]=="activation") ) {
      nBaseParam <- 3
    } else {
      nBaseParam <- 4
    }
  }

  nParamk <- nBaseParam + ngParamk
  funcName <- paste(model$Transform, "Transform", sep="")
  func <- get(funcName, mode="function")

  for ( i in seq(length=model$numGenes) ) {
    if ( any(grep("bTransform", names(model))) ) {
      if ( is.na(model$bTransform) ) {
        model$B[i] <- params[endVal+nParamk*(i-1)+1]
      } else {
        model$B[i] <- func(params[endVal+nParamk*(i-1)+1], "atox")
      }
    } else {
      model$B[i] <- func(params[endVal+nParamk*(i-1)+1], "atox")
    }
    model$S[i] <- func(params[endVal+nParamk*(i-1)+2], "atox")
    model$D[i] <- func(params[endVal+nParamk*(i-1)+3], "atox")

    if ( nBaseParam == 4 )
      if ( any(grep("alphaTransform", names(model))) ) {
        if ( is.na(model$alphaTransform) ) {
          model$alpha[i] <- params[endVal+nParamk*(i-1)+4]   
        } else {
          model$alpha[i] <- func(params[endVal+nParamk*(i-1)+4], "atox")
        }
      } else {
        model$alpha[i] <- func(params[endVal+nParamk*(i-1)+4], "atox")
      }

    if ( model$ngParam>0 ) {
      startV <- endVal+nParamk*(i-1)+nBaseParam+1
      endV <- endVal+nParamk*(i-1)+nParamk
      if ( length(model$gParam) > 1 ) {
        model$gParam[,i] <- func(params[startV:endV], "atox")
      } else if ( length(model$gParam) == 1 ) {
        model$gParam <- func(params[startV:endV], "atox")
      }
    }
  }

  if ( model$includeNoise ) {
    noiseSd <- func(params[(length(params)-model$numGenes+1):length(params)], "atox")
    model$noiseVar <- noiseSd*noiseSd
  }
  
  model <- gpsimMapUpdateKernels(model)
  model <- gpsimMapUpdatePosteriorCovariance(model)
  model$updateW <- TRUE
  
  return (model)
}



gpsimMapObjective <- function (param, model) {
  if ( any(grep("comp", names(model))) ) {
    Nrep <- length(model$comp)
    ll <- array(0, Nrep)
    for ( rep in seq(length=Nrep) ) {
      options <- optimiFdefaultOptions()
      model$comp[[rep]] <- gpsimMapExpandParam(model$comp[[rep]], param)
      model$comp[[rep]] <- gpsimMapUpdateF(model$comp[[rep]], options)
      ll[rep] <- gpsimMapLogLikelihood(model$comp[[rep]])
    }
  }

  f <- -sum(ll)    

  return (f)
}


optimiFdefaultOptions <- function () {
  options <- list(maxit=100, tolf=1e-4, tol=1e-4, display=FALSE)
  return (options)
}



gpsimMapUpdateF <- function (model, options) {
  if ( !any(grep("updateF",names(model))) ) {
    updateF <- TRUE
  } else if ( model$updateF ) {
    updateF <- TRUE
  } else {
    updateF <- FALSE
  }
  
  if ( updateF ) {
    cat("\n Updating F ...... \n")
    display <- options$display
    tolf <- options$tolf
    tol <- options$tol
    iters <- options$maxit
    
    f <- gpsimMapFunctionalExtractParam(model)
    
    if ( any(grep("proteinPrior",names(model))) ) {
      model$consLambda <- 5
      maxLambda <- 8
      nCons <- length(model$proteinPrior$values)
    } else {
      model$consLambda <- 0
      maxLambda <- 1
    }
    
    llold <- gpsimMapFunctionalLogLikelihood(model)
    
    for ( j in seq(length=maxLambda) ) {
      ## cat("Lambda=", model$consLambda, "\n")
      flag <- 0
      for ( i in seq(length=iters) ) {
        ll <- llold-1
        gf <- gpsimMapFunctionalLogLikeGradients(model)
        if ( model$consLambda ) {
          consMat <- matrix(0, dim(model$invCovf)[1], dim(model$invCovf)[2])
          for ( k in seq(length=nCons) ) {
            ftimeIndex <- grep(TRUE, (model$proteinPrior$times[k]-model$mapt)==0)
            consMat[ftimeIndex, ftimeIndex] <- model$consLambda          
          }
          
          # hf <- solve(model$invCovf+consMat)
          hf <- jitCholInv(model$invCovf+consMat)$invM
        } else {
          hf <- model$covf
        }
        
        newDrn <- t(hf%*%gf)[1,]
        factor <- 1
        count <- 0
        fold <- f
        
        while ( ll<llold ) {
          f <- fold+factor*newDrn
          model <- gpsimMapFunctionalExpandParam(model, f)
          ll <- gpsimMapFunctionalLogLikelihood(model)
          lldiff <- ll-llold
          count <- count+1
          if ( count>0 & display )
            cat("gpsimMapUpdateF, lldiff: ", round(lldiff, digits=4), ", factor ", factor, "\n")
          
          factor <- factor/2
          if ( lldiff<tol & max(abs(fold-f))<tolf ) {
            flag <- 1
            break
          }
        }
        
        if ( display )
          cat("Iteration ", i, ", log likelihood ", ll, "\n")
        
        llold <- ll
        if (flag == 1)
          break
      }
      model$consLambda <- model$consLambda*4
    }
    
    if ( all(model$f==0) )
      warning("TF is zero!")
  } else {
    if ( options$display )
      cat("\n F is not updated! \n")
  }
  return (model)
}



gpsimMapFunctionalExtractParam <- function (model) {
  return (model$f)
}



gpsimMapFunctionalLogLikelihood <- function (model) {
  ll <- gpsimMapLogLikelihood(model) - 0.5*model$logDetCovf

  if ( any(names(model)=="proteinPrior") ) {
    nCons <- length(model$proteinPrior$values)
    for ( k in seq(length=nCons) ) {
      ftimeIndex <- grep(TRUE, (model$proteinPrior$times[k]-model$mapt)==0)
      ll <- ll-0.5*model$consLambda*(model$f[ftimeIndex]-model$proteinPrior$values[k])^2
    }
  } 

  return (ll)
}



gpsimMapLogLikelihood <- function (model) {
  numData <- length(model$t)
  ll <- t(model$f) %*% model$invK %*% model$f - model$logDetCovf + model$logDetK

  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
  } else {
    yvar <- model$yvar
  }

  for ( i in seq(length=numData) )
    for ( j in seq(length=model$numGenes) ) {
      ind <- i + (j-1)*numData
      betaij <- 1/yvar[ind]
      factor <- (model$ypred[model$timesIndex[i], j] - model$y[ind])
      ll <- ll + factor*factor*betaij - log(betaij)  
    }

  ll <- ll + numData*model$numGenes*log(2*pi)
  ll <- -0.5*ll
  
  return (ll)
}



gpsimMapFunctionalLogLikeGradients <- function (model) {

  gdata <- array(0, model$numMapPts)
  numData <- length(model$t)

  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
  } else {
    yvar <- model$yvar
  }

  for ( k in seq(model$timesIndex[1]+1, model$numMapPts) ) {
    temp <- 0
    for ( i in seq(length=numData) ) {
      arg <- model$t[i]-model$mapt[k]
      if ( arg>=0 )
        for ( j in seq(length=model$numGenes) ) {
          if ( model$ngParam ) {
            gInd <- j
          } else {
            gInd <- 1
          }
          ind <- (i + (j-1)*numData)
          betaij <- 1/yvar[ind]
          factor <- (model$ypred[model$timesIndex[i], j]-model$y[ind])*betaij
          temp <- temp+factor*as.matrix(model$gGrad)[k,gInd]*exp(-model$D[j]*arg)*model$S[j]
        }
    }
    gdata[k] <- -temp*model$step
  }

  g <- gdata - t(model$invK%*%model$f)[1,]

  gCons <- array(0, length(g))

  if ( any(grep("proteinPrior",names(options))) ) {
    nCons <- length(model$proteinPrior$values)
    for ( i in seq(length=nCons) ) {
      ftimeIndex <- grep(TRUE, (model$proteinPrior$times[i]-model$mapt)==0)
      g[ftimeIndex] <- g[ftimeIndex]-model$consLambda*(model$f[ftimeIndex]-model$proteinPrior$values[i])
    }
  }   

  return (g)
}



gpsimMapFunctionalExpandParam <- function (model, f) {
  model$f <- f
  model <- gpsimMapUpdateG(model)
  model <- gpsimMapUpdateYpred(model)
  model <- gpsimMapFunctionalUpdateW(model)
  model <- gpsimMapUpdatePosteriorCovariance(model)
  return (model)
}


  
gpsimMapGradients <- function (param, model) {

  g <- array(0, length(param))
  Nrep <- length(model$comp)

  dg <- list()
   
  for ( rep in seq(length=Nrep) ) {
    options <- optimiFdefaultOptions()
    model$comp[[rep]] <- gpsimMapExpandParam(model$comp[[rep]], param)
    model$comp[[rep]] <- gpsimMapUpdateF(model$comp[[rep]], options)
    dg[[rep]] <- gpsimMapLogLikeGradients(model$comp[[rep]])
    g <- g - dg[[rep]]
  }

  return (g)
}



gpsimMapLogLikeGradients <- function (model) {

  covGrad <- model$invK %*% (model$f%*%t(model$f)+model$covf) %*% model$invK - model$invK
  covGrad <- covGrad*0.5

  gkern <- kernGradient(model$kern, model$mapt, covGrad)
  gmodel <- gpsimMapMarginalLikeliGradient(model)

  g <- c(gkern, gmodel)

  if ( all(model$nonLinearity=="linear") ) {
    g1 <- 0
  } else {
    g1 <- gpsimMapLikeliGradientImplicit(model)
  }

  g <- g + g1

  if ( any(grep("fix",names(model))) ) 
    g[model$fix$index] <- 0

  return (g)
}



gpsimMapMarginalLikeliGradient <- function (model) {

  numData <- length(model$t)
  ParamL <- length(model$B)

  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
    dlogPdn <- matrix(0, ParamL, 1)
    gn <- array(0, model$numGenes)
  } else {
    yvar <- model$yvar
  }

  if ( any(grep("includeRepression", names(model))) )
    if ( model$includeRepression ) {
      galpha <- array(0, model$numGenes)
      dlogPdalpha <- matrix(0, ParamL, 1)
    }

  gB <- array(0, model$numGenes)
  gD <- array(0, model$numGenes)
  gS <- array(0, model$numGenes)

  dlogPdB <- matrix(0, ParamL, 1)
  dlogPdD <- matrix(0, ParamL, 1)
  dlogPdS <- matrix(0, ParamL, 1)

  if ( model$ngParam > 0 ) {
    ngParamk <- model$ngParam/model$numGenes
    dlogPdgParam <- matrix(0, ParamL, ngParamk)
    ggParam <- matrix(0, ngParamk, ParamL)
  }

  invC <- model$covf

  for ( k in seq(length=model$numGenes) ) {
    for ( i in seq(length=numData) ) {
      ind <- i + (k-1)*numData
      betaik <- 1/yvar[ind]
      factor <- (model$ypred[model$timesIndex[i], k]-model$y[ind])*betaik

      dx <- gpsimXGradient(model, i, k)
      dxdB <- dx$dxdB
      dxdD <- dx$dxdD
      dxdS <- dx$dxdS
      dxdalpha <- dx$dxdalpha
      dxdgParam <- dx$dxdgParam

      dlogPdB[k] <- dlogPdB[k]-factor*dxdB

      if ( any(grep("includeRepression", names(model))) )
        if ( model$includeRepression )
          dlogPdalpha[k] <- dlogPdalpha[k]-factor*dxdalpha
                    
      dlogPdD[k] <- dlogPdD[k]-factor*dxdD
      dlogPdS[k] <- dlogPdS[k]-factor*dxdS

      if ( model$ngParam > 0 )
        dlogPdgParam[k,] <- dlogPdgParam[k,]-factor*dxdgParam

      if ( model$includeNoise )
        dlogPdn[k] <- dlogPdn[k]+sqrt(model$noiseVar[k])*factor*factor-sqrt(model$noiseVar[k])*betaik
    }

    dW <- gpsimMapWGradient(model, k)
    dWdB <- dW$dWdB
    dWdD <- dW$dWdD
    dWdS <- dW$dWdS
    dWdalpha <- dW$dWdalpha
    dWdgParam <- dW$dWdgParam
    dWdn <- dW$dWdn
 
    gB[k] <- -0.5*matrixTrace(invC%*%dWdB)+dlogPdB[k]
    gD[k] <- -0.5*matrixTrace(invC%*%dWdD)+dlogPdD[k]
    gS[k] <- -0.5*matrixTrace(invC%*%dWdS)+dlogPdS[k]

    if ( any(grep("includeRepression", names(model))) )
      if ( model$includeRepression )
        galpha[k] <- -0.5*matrixTrace(invC%*%dWdalpha)+dlogPdalpha[k]

    if ( model$ngParam > 0 )
      for ( gParamIndex in seq(length=ngParamk) )
        ggParam[gParamIndex,k] <- -0.5*matrixTrace(invC%*%dWdgParam[[gParamIndex]])+dlogPdgParam[k,gParamIndex]
    if ( model$includeNoise )
      gn[k] <- -0.5*matrixTrace(invC%*%dWdn) + dlogPdn[k]
  }

  funcName <- paste(model$Transform, "Transform", sep="")
  func <- get(funcName, mode="function")
  g <- NULL

  for ( i in seq(length=model$numGenes) ) {
    if ( any(grep("bTransform", names(model))) ) {
      if ( is.na(model$bTransform) ) {
        g <- c(g, gB[i])
      } else {
        g <- c(g, gB[i]*func(model$B[i], "gradfact"))
      }
    } else {
      g <- c(g, gB[i]*func(model$B[i], "gradfact"))
    }

    g <- c(g, gS[i]*func(model$S[i], "gradfact"), gD[i]*func(model$D[i], "gradfact"))

    if ( any(grep("includeRepression", names(model))) )
      if ( model$includeRepression )
        if ( any(grep("alphaTransform", names(model))) ) {
          if ( is.na(model$alphaTransform) ) {
            g <- c(g, galpha[i])
          } else {
            g <- c(g, galpha[i]*func(model$alpha[i], "gradfact"))
          }
        } else {
          g <- c(g, galpha[i]*func(model$alpha[i], "gradfact"))
        }

    if ( model$ngParam>0 )
      g <- c(g, ggParam[,i]*func(as.matrix(model$gParam)[,i], "gradfact"))
  }

  if ( model$includeNoise )
    g <- c(g, gn*func(sqrt(model$noiseVar), "gradfact"))

  return (g)
}



gpsimMapLikeliGradientImplicit <- function (model) {
  g1 <- array(0, model$kern$nParams)
  g2 <- array(0, model$numParams-model$kern$nParams)

  covGrad <- model$covf%*%model$invK
  funcGrad <- model$invK%*%model$f

  kernGrad <- matrix(0, model$numMapPts, model$kern$nParams)
  dfuncGrad <- matrix(0, model$numMapPts, model$numParams-model$kern$nParams)

  for ( ftimeIndex in seq(length=model$numMapPts) ) {
    kernGrad[ftimeIndex,] <- kernGradient(model$kern, model$mapt, covGrad)*funcGrad[ftimeIndex]
    dfuncGrad[ftimeIndex,] <- gpsimMapFunctionalLikeGrad2(model, ftimeIndex)
  }

  modelParamGrad <- covGrad%*%model$K%*%dfuncGrad

  dWdf <- list()
  dlogPdf <- array()

  for ( ftimeIndex in seq(length=model$numMapPts) ) {
    dWdf[[ftimeIndex]] <- gpsimMapFunctionalWGradient(model, ftimeIndex)
    dlogPdf[ftimeIndex] <- -0.5*matrixTrace(covGrad%*%model$K%*%dWdf[[ftimeIndex]])
    g1 <- g1 + dlogPdf[ftimeIndex]*kernGrad[ftimeIndex,]
    g2 <- g2 + dlogPdf[ftimeIndex]*modelParamGrad[ftimeIndex,]
  }

  ## Check if model parameters are being optimised in a transformed space
  params <- gpsimMapExtractParam(model, only.values=FALSE)
  paramvec <- params$values
  names <- params$names

  if ( any(grep("bTransform", names(model))) )
    if ( is.na(model$bTransform) ) { 
      Bindex <- grep("Basal", names)
      paramvec[Bindex] <- 0
    }

  if ( any(grep("alphaTransform", names(model))) )
    if ( is.na(model$alphaTransform) ) { 
      alphaIndex <- grep("Param", names)
      paramvec[alphaIndex] <- 0
    }   

  modelParams <- paramvec[(model$kern$nParams+1):length(paramvec)]

  if ( !is.na(model$Transform) ) {
    funcName <- paste(model$Transform, "Transform", sep="")
    func <- get(funcName, mode="function")
    modelParams <- func(modelParams, "atox")
    g2 <- g2*func(modelParams, "gradfact")
  }

  g <- c(g1, g2)  

  return (g)
}



gpsimXGradient <- function (model, i, j, option=1) {

  dxdB <- 1/model$D[j]
  endPoint <- model$timesIndex[i]

  dxdD <- 0
  dxdS <- 0
  dxdalpha <- array()
  dxdgParam <- matrix()

  if ( model$ngParam > 0 ) {
    ngParamk <- model$ngParam/model$numGenes
    dxdgParam <- matrix(0, 1, ngParamk)
    gInd <- j
  } else {
    gInd <- 1
  }

  if ( option == 1 ) {
    for ( m in seq(length=endPoint) ) {
      arg <- model$t[i]-model$mapt[m]
      if ( arg>=0 ) {
        dxdD <- dxdD + as.matrix(model$g)[m,gInd]*arg*exp(-model$D[j]*arg+log(model$step)+log(model$S[j]))
        dxdS <- dxdS + exp(-model$D[j]*arg+log(model$step))*as.matrix(model$g)[m,gInd]
        if ( model$ngParam > 0 ) 
          for ( gParamInd in seq(length=ngParamk) ) 
            dxdgParam[gParamInd] <- dxdgParam[gParamInd] + exp(-model$D[j]*arg+log(model$step)+log(model$S[j]))*as.matrix(model$dg)[m,gInd]
      }
    }

  } else if ( option==2 ) {  
    integral <- exp(model$D[j]*model$mapt[1:endPoint])
    if ( endPoint%%2 ) {         
      integral1 <- integral
      integral1[2*seq(length=floor((length(integral1)-1)/2))] <- 4*integral[2*seq(length=floor((length(integral1)-1)/2))]
      if ( length(integral1)>1 )
        integral1[2*seq(length=floor((length(integral1)-2)/2))+1] <- 2*integral[2*seq(length=floor((length(integral1)-2)/2))+1]
    } else { 
      integral1 <- integral
      integral1[1] <- 3*integral[1]
      integral1[2*seq(length=floor((length(integral1)-1)/2))+1] <- 4*integral[2*seq(length=floor((length(integral1)-1)/2))]
      if ( length(integral1)>1 )
        integral1[2*seq(length=floor((length(integral1)-2)/2))+2] <- 2*integral[2*seq(length=floor((length(integral1)-2)/2))+1]
    }
    
    dxdD <- exp(log(model$S[j])-model$D[j]*model$t[i]+log(model$step))/3*sum(model$g[1:endPoint,gInd]*integral1[1:endPoint]*(model$t[i]-model$mapt[1:endPoint]))
  
    dxdS <- exp(-model$D[j]*model$t[i]+log(model$step))/3*sum(as.matrix(model$g)[1:endPoint,gInd]*integral1[1:endPoint])
    
    if ( model$ngParam > 0 )
      for ( gParamInd in seq(length=ngParamk) )
        dxdgParam[gParamInd] <- exp(log(model$S[j])-model$D[j]*model$t[i]+log(model$step))/3*sum(as.matrix(model$dg)[1:endPoint, gInd]*integral1[1:endPoint])
  }
  
  dxdD <- -model$B[j]/(model$D[j]*model$D[j])-dxdD

  if ( any(grep("isGroupNonlinearity", names(model))) )
    if ( model$nonLinearity[j]=="repression" ) {
      dxdalpha <- exp(-model$D[j]*model$t[i])
      dxdD <- -model$t[i]*model$alpha[j]*exp(-model$D[j]*model$t[i])+dxdD
    }

  dx <- list(dxdB=dxdB, dxdD=dxdD, dxdS=dxdS, dxdalpha=dxdalpha, dxdgParam=dxdgParam)  

  return (dx)
}



gpsimMapWGradient <- function(model, k) {

  intPoints <- (model$timesIndex[1]+1):model$numMapPts
  step2 <- model$step*model$step
  S2 <- model$S*model$S
  numData <- length(model$t)

  w1 <- dim(model$W)[1]
  w2 <- dim(model$W)[2]

  dWdB <- matrix(0, w1, w2)
  dWdD <- matrix(0, w1, w2)
  dWdS <- matrix(0, w1, w2)
  dWdalpha <- matrix()
  dWdgParam <- matrix()
  dWdn <- matrix()

  if ( model$ngParam > 0 ) {
    ngParamk <- model$ngParam/model$numGenes
    dWdgParam <- list()
    for ( i in seq(length=ngParamk) )
      dWdgParam[[i]] <- matrix(0, w1, w2)
    gInd <- k
  } else {
    gInd <- 1
  }

  if ( any(grep("isGroupNonlinearity", names(model))) )
    if ( model$nonLinearity[k]=="repression" )
      dWdalpha <- matrix(0, w1, w2)

  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
    dWdn <- matrix(0, w1, w2)
  } else {
    yvar <- model$yvar
  }

  for ( p in intPoints )
    for ( i in seq(length=numData) ) {
      arg <- model$t[i]-model$mapt[p]
      if ( arg >= 0 ) {
        ind <- i + (k-1)*numData
        betaik <- 1/yvar[ind]

        dx <- gpsimXGradient(model, i, k)
        dxdB <- dx$dxdB
        dxdD <- dx$dxdD
        dxdS <- dx$dxdS
        dxdalpha <- dx$dxdalpha
        dxdgParam <- dx$dxdgParam
        
        dWdB[p, p] <- dWdB[p, p]+betaik*as.matrix(model$gGrad2)[p,gInd]*dxdB*exp(-model$D[k]*arg+log(model$S[k]) +log(model$step))

        if ( any(grep("isGroupNonlinearity", names(model))) )
          if ( model$nonLinearity[k]=="repression" )
            dWdalpha[p,p] <- dWdalpha[p, p]+betaik*as.matrix(model$gGrad2)[p,gInd]*dxdalpha*exp(-model$D[k]*arg+log(model$S[k]) +log(model$step))
                 
        factor <- model$ypred[model$timesIndex[i], k]-model$y[ind]

        dWdD[p, p] <- dWdD[p, p]+model$step*betaik*as.matrix(model$gGrad2)[p,gInd]*(dxdD-factor*arg)*exp(-model$D[k]*arg+log(model$S[k]))

        dWdS[p, p] <- dWdS[p, p]+model$step*betaik*exp(-model$D[k]*arg)*(dxdS*model$S[k]*as.matrix(model$gGrad2)[p,gInd]+factor*as.matrix(model$gGrad2)[p,gInd])
        
        if ( model$ngParam > 0 )
          for ( gParamInd in seq(length=ngParamk) )
            dWdgParam[[gParamInd]][p, p] <- dWdgParam[[gParamInd]][p, p]+model$step*betaik*exp(-model$D[k]*arg)*(dxdgParam[gParamInd]*model$S[k]*model$gGrad2[p,gInd]+factor*model$S[k]*model$dggrad2[p,gInd])
        
        if ( model$includeNoise )
          dWdn[p, p] <- dWdn[p, p]-2*sqrt(model$noiseVar[k])*betaik^2*factor*as.matrix(model$gGrad2)[p,gInd]*exp(-model$D[k]*arg+log(model$step)+log(model$S[k]))
      }
    }

  for ( p in intPoints )
    for ( q in intPoints )
      for ( i in seq(length=numData) ) {
        arg1 <- model$t[i]-model$mapt[p]
        arg2 <- model$t[i]-model$mapt[q]
        if ( arg1 >= 0 & arg2 >= 0 ) {
          ind <- i + (k-1)*numData
          betaik <- 1/yvar[ind]

          dWdD[p,q] <- dWdD[p,q]-betaik*as.matrix(model$gGrad)[p,gInd]*as.matrix(model$gGrad)[q,gInd]*(arg1+arg2)*exp(-model$D[k]*(arg1+arg2)+log(S2[k])+log(step2))

          dWdS[p,q] <- dWdS[p,q]+2*betaik*model$S[k]*as.matrix(model$gGrad)[q,gInd]*as.matrix(model$gGrad)[p,gInd]*exp(-model$D[k]*(arg1+arg2)+log(step2))

          if ( model$ngParam > 0 )
            for ( gParamInd in seq(length=ngParamk) )
              dWdgParam[[gParamInd]][p, q] <- dWdgParam[[gParamInd]][p, q]+betaik*(model$dggrad[q,gInd]*model$gGrad[p,gInd]+model$gGrad[q,gInd]*model$dggrad[p,gInd])*exp(-model$D[k]*(arg1+arg2)+log(S2[k])+log(step2))

          if ( model$includeNoise )
            dWdn[p, q] <- dWdn[p, q]-2*sqrt(model$noiseVar[k])*betaik^2*as.matrix(model$gGrad)[p,gInd]*as.matrix(model$gGrad)[q,gInd]*exp(-model$D[k]*(arg1+arg2)+log(step2)+log(S2[k]))
        }
      }

  dW <- list(dWdB=dWdB, dWdD=dWdD, dWdS=dWdS, dWdalpha=dWdalpha, dWdgParam=dWdgParam, dWdn=dWdn)  

  return (dW)
}
                             


gpsimMapFunctionalWGradient <- function (model, l, option=1) {

  intPoints <- (model$timesIndex[1]+1):model$numMapPts
  step2 <- model$step*model$step
  S2 <- model$S*model$S
  numData <- length(model$t)

  dWdf <- matrix(0, dim(model$W)[1], dim(model$W)[2])

  if ( model$includeNoise ) {
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
  } else {
    yvar <- model$yvar
  }
  
  for ( p in intPoints ) {
    ## part1: p=l or q=l
    for ( q in intPoints )
      for ( i in seq(length=numData) ) {
        arg1 <- model$t[i]-model$mapt[p]
        arg2 <- model$t[i]-model$mapt[q]
        if ( arg1 >= 0 & arg2 >= 0 ) {
          if ( p==l )
            for ( j in seq(length=model$numGenes) ) {
              if ( model$ngParam ) {
                gInd <- j
              } else {
                gInd <- 1
              }

              ind <- i + (j-1)*numData
              betaij <- 1/yvar[ind]
              dWdf[p,q] <- dWdf[p,q]+betaij*model$gGrad2[p,gInd]*model$gGrad[q,gInd]*exp(-model$D[j]*(arg1+arg2)+log(S2[j])+log(step2))
            }
             
          if ( q==l )
            for ( j in seq(length=model$numGenes) ) {
              if ( model$ngParam ) {
                gInd <- j
              } else {
                gInd <- 1
              }
              ind <- i + (j-1)*numData
              betaij <- 1/yvar[ind]
              dWdf[p,q] <- dWdf[p,q]+betaij*model$gGrad2[q,gInd]*model$gGrad[p,gInd]*exp(-model$D[j]*(arg1+arg2)+log(S2[j])+log(step2))
            }
        }
      }

    ## part2: p=q
    for ( i in seq(length=numData) ) {
      arg1 <- model$t[i]-model$mapt[p]

      if ( option==1 ) {
        arg2 <- model$t[i]-model$mapt[l]
        if ( arg1 >= 0 & arg2 >= 0 ) {  
          for ( j in seq(length=model$numGenes) ) {
            if ( model$ngParam ) {
              gInd <- j
            } else {
              gInd <- 1
            }
            
            ind <- i + (j-1)*numData
            betaij <- 1/yvar[ind]
            dWdf[p,p] <- dWdf[p,p]+betaij*model$gGrad[l,gInd]*model$gGrad2[p,gInd]*exp(-model$D[j]*(arg1+arg2)+log(S2[j])+log(step2))            
          }
        }     
      } else if ( option==2 ) {
        if ( arg1 >= 0 ) {  
          for ( j in seq(length=model$numGenes) ) {
            if ( model$ngParam ) {
              gInd <- j
            } else {
              gInd <- 1
            }
            ind <- i + (j-1)*numData
            betaij <- 1/yvar[ind]
            
            dxdf <- gpsimMapFunctionalXGradient(model, i, j, l) 
            dWdf[p,p] <- dWdf[p,p]+betaij*dxdf*model$gGrad2[p,gInd]*exp(-model$D[j]*arg1+log(S2[j])+log(step2))
          }
        }
      }
    }
  }
  
  ## part3: p=q=l
  temp <- 0
  for ( i in seq(length=numData) ) {
    arg <- model$t[i]-model$mapt[l]
    if ( arg >= 0  & l > 1 )
      for ( j in seq(length=model$numGenes) ) {
        if ( model$ngParam ) {
          gInd <- j
        } else {
          gInd <- 1
        }
         
        ind <- i + (j-1)*numData
        betaij <- 1/yvar[ind]
        factor <- (model$ypred[model$timesIndex[i], j]-model$y[ind])*betaij
        temp <- temp+factor*model$gGrad3[l,gInd]*exp(-model$D[j]*arg+log(model$S[j])+log(model$step))
      }
  }

  dWdf[l, l] <- dWdf[l, l] + temp  

  return (dWdf)
}



gpsimMapFunctionalLikeGrad2 <- function (model, p) {

  dfuncGrad <- array(0, model$numParams-model$kern$nParams)
  gB <- array(0, model$numGenes)
  gD <- array(0, model$numGenes)
  gS <- array(0, model$numGenes)
  ngParamk <- model$ngParam/model$numGenes
  numData <- length(model$t)

  if ( model$includeNoise ) {
    gn <- array(0, model$numGenes)
    noiseMat <- matrix(rep(model$noiseVar, each=numData), numData, model$numGenes)
    yvar <- model$yvar + noiseMat
  } else {
    yvar <- model$yvar
  }

  if ( any(grep("includeRepression", names(model))) )
    if ( model$includeRepression ) {
      nBaseParam <- 4
      galpha <- array(0, model$numGenes)
    } else {
      nBaseParam <- 3
    } else {
      nBaseParam <- 3
    }
  
  if ( model$ngParam > 0 )
    ggParam <- matrix(0, model$numGenes, ngParamk)
  
  for ( k in seq(model$numGenes) ) {
    if ( model$ngParam ) {
      gInd <- k
    } else {
      gInd <- 1
    }
  
    for ( i in seq(length=numData) ) {
      arg <- model$t[i]-model$mapt[p]
      if ( arg >= 0 ) {
        ind <- i + (k-1)*numData
        betaik <- 1/yvar[ind]

        dx <- gpsimXGradient(model, i, k)
        dxdB <- dx$dxdB
        dxdD <- dx$dxdD
        dxdS <- dx$dxdS
        dxdalpha <- dx$dxdalpha
        dxdgParam <- dx$dxdgParam
        
        gB[k] <- gB[k] - betaik*model$gGrad[p,gInd]*dxdB*exp(-model$D[k]*arg+log(model$S[k])+log(model$step))

        factor <- model$ypred[model$timesIndex[i], k] - model$y[ind]
        gD[k] <- gD[k] - betaik*model$gGrad[p,gInd]*(dxdD-arg*factor)*exp(-model$D[k]*arg+log(model$S[k])+log(model$step))
        
        gS[k] <- gS[k] - betaik*model$gGrad[p,gInd]*(dxdS*model$S[k]+factor)*exp(-model$D[k]*arg+log(model$step))

        if ( any(grep("includeRepression", names(model))) )
          if ( model$includeRepression )
            galpha[k] <- galpha[k] - betaik*model$gGrad[p,gInd]*dxdalpha*exp(-model$D[k]*arg+log(model$S[k])+log(model$step))
             
        if ( model$ngParam > 0 )
          ggParam[k,] <- ggParam[k,] - betaik*(model$gGrad[p,gInd]*dxdgParam+factor*model$dggrad[p, gInd])*exp(-model$D[k]*arg+log(model$S[k])+log(model$step))
        
        if ( model$includeNoise )
          gn[k] <- gn[k] + 2*sqrt(model$noiseVar[k])*betaik^2*model$gGrad[p,gInd]*model$S[k]*factor*exp(-model$D[k]*arg+log(model$step))
      }
    }
  }
  
  geneParam <- nBaseParam + ngParamk  
  startPoint <- 1
  endPoint <- geneParam
  
  for ( k in seq(length=model$numGenes) ) {
    dfuncGrad[startPoint:(startPoint+2)] <- c(gB[k], gS[k], gD[k])
  
    if ( any(grep("includeRepression", names(model))) )
      if ( model$includeRepression )
        dfuncGrad[startPoint+3] <- galpha[k]
   
    if ( model$ngParam > 0 )
      dfuncGrad[(startPoint+nBaseParam):endPoint] <- ggParam[k,]
        
    startPoint <- startPoint + geneParam
    endPoint <- endPoint + geneParam
  }

  if ( model$includeNoise )
    dfuncGrad[(length(dfuncGrad)-length(gn)+1):length(dfuncGrad)] = gn
  
  return (dfuncGrad)
}



gpsimMapFunctionalXGradient <- function (model, i, j, l) {

  endPoint <- model$timesIndex[i]
  integral <- array(0, model$numMapPts)
  
  if ( model$ngParam > 0 ) {
    gInd <- j
  } else {
    gInd <- 1
  }

  integral[1:endPoint] <- exp(model$D[j]*model$mapt[1:endPoint])
  if ( endPoint%%2 ) {         
    integral1 <- integral
    integral1[2*seq(length=floor((length(integral1)-1)/2))] <- 4*integral[2*seq(length=floor((length(integral1)-1)/2))]
    if ( length(integral1)>1 )
        integral1[2*seq(length=floor((length(integral1)-2)/2))+1] <- 2*integral[2*seq(length=floor((length(integral1)-2)/2))+1]
  } else { 
    integral1 <- integral
    integral1[1] <- 3*integral[1]
    integral1[2*seq(length=floor((length(integral1)-1)/2))+1] <- 4*integral[2*seq(length=floor((length(integral1)-1)/2))]
    if ( length(integral1)>1 )
      integral1[2*seq(length=floor((length(integral1)-2)/2))+2] <- 2*integral[2*seq(length=floor((length(integral1)-2)/2))+1]
  }

  dxdf <- exp(log(model$S[j])-model$D[j]*model$t[i]+log(model$step))/3*model$gGrad[l,gInd]*integral1[l]
  
  return (dxdf)
}



gpsimMapUpdateYpredVar <- function (model) {

  iter <- 100
  predModel <- list()
  ypred <- list()
  for ( j in seq(length=model$numGenes) ) {
    ypred[[j]] <- matrix(, length(model$f), iter)
  }
  for ( i in 1:iter ) {
    f <- model$f + sqrt(model$varf)*rnorm(length(model$f))
    predModel$samp[[i]] <- gpsimMapFunctionalExpandParam(model, f)
    
    for ( j in seq(length=model$numGenes) ) {
       ypred[[j]][,i] <- predModel$samp[[i]]$ypred[,j]
    }
  }

  ypredVar <- matrix(, length(model$f), model$numGenes)
  for ( j in seq(length=model$numGenes) ) {
    if ( model$includeNoise ) {
      ypredVar[,j] <- apply(ypred[[j]], 1, var) + model$noiseVar[j]
    } else {
      ypredVar[,j] <- apply(ypred[[j]], 1, var)
    }    
  }

  model$ypredVar <- ypredVar
  return (model)
}



cgpsimMapOptimise <- function (model, options, ...) {
  if ( any(grep("optimiser", names(options))) ) {
    funcName <- paste(options$optimiser, "optim", sep="")
  } else {
    funcName <- "CGoptim"
  }
  func <- get(funcName, mode="function")

  params <- modelExtractParam(model$comp[[1]])
  newParams <- func(params, fn=gpsimMapObjective, grad=gpsimMapGradients, options, model)

  Nrep <- length(model$comp)
  for ( i in seq(length=Nrep) ) {
    optimiFoptions <- optimiFdefaultOptions()
    model$comp[[i]] <- gpsimMapExpandParam(model$comp[[i]], newParams$xmin)
    model$comp[[i]] <- gpsimMapUpdateF(model$comp[[i]], optimiFoptions)
    model$comp[[i]] <- gpsimMapUpdateYpredVar(model$comp[[i]])
  }
  
  model$llscore <- newParams$objective
  
  if ( funcName == "CGoptim" )
    model$lnSchFail <- newParams$lnSchFail
 
  return (model)
}



gpsimMapBarencoResults <- function (model, expType, expNo, scale=1, option=1) {
  geneNames <- c("DDB2", "hPA26", "TNFRSF10b", "p21", "BIK")
  order <- c(1,5,3,4,2)

  modelB <- model$comp[[1]]$B*scale
  if ( model$comp[[1]]$ngParam ) {
    modelS <- (model$comp[[1]]$S/(1+model$comp[[1]]$gParam)[1,])*scale
  } else {
    modelS <- model$comp[[1]]$S*scale
  }
  modelS <- modelS/modelS[4]

  if ( model$comp[[1]]$ngParam ) {
    nonLinearity <- "multi"
  } else {
    nonLinearity <- model$comp[[1]]$nonLinearity
  }

  if ( option == 1 ) {

  for ( j in seq(along=model$comp) ) {
  
    if ( model$comp[[1]]$ngParam ) {
      scalePred <- sqrt(var(exp(model$comp[[j]]$f)))
    } else {
      scalePred <- sqrt(var(model$comp[[j]]$f)) 
    }
  
    ## Info from Martino Paper:
    ## 'True f' from Figure 3.
    truef <- c(0, 1.6, 2.6, 2.5, 2.6, 1.6, 0.9)
    truef <- truef/sqrt(var(truef))*scalePred
    
    ## Figure 2(a) histograms;
    B <- c(2.6, 1.5, 0.5, 0.2, 1.35) ## From Martino paper ... but don't know the scale
    B <- B/mean(B)*mean(modelB)      ## do a rough rescaling so that the scales match.
    S <- c(3, 0.8, 0.7, 1.8, 0.7)/1.8  
    D <- c(1.2, 1.6, 1.75, 3.2, 2.3)*0.8/3.2; 
  
    ## Martino f from Figure 2(b), again measured with a ruler.
    barencof <- t(matrix(c(0.0000000, 200.5201100, 355.5216125, 205.7574913, 135.0911372, 
              145.1080997, 130.7046969, 0.0000000, 184.0994134, 308.4759200, 
              232.1775328, 153.6595161, 85.7272235,  168.0910562, 0.0000000, 
              230.2262511, 337.5994811, 276.9416540, 164.5044287, 127.8653452, 173.6112139), 7, 3))
    barencof <- barencof/(1.8*mean(S))*mean(modelS)
    barencof <- barencof/(sqrt(apply(barencof, 1, var))%*%matrix(1,1,7))*scalePred

    x11()
    if ( nonLinearity == "linear" ) {
      fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
      fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
      plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
    } else if ( nonLinearity == "multi" ) {
      fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
      fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
      plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
    } else {
      funcName <- model$comp[[j]]$nonLinearity
      func <- get(funcName, mode="function")
      fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
      fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
      plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
    }

    for ( index in seq(length=model$comp[[j]]$numGenes) ) {
      x11()
      xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]))+0.1
      xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]))-0.1
      plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
      title(paste("Replica ", j, ": ", index, "-th Gene", sep=""))
      lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
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
 
  for ( j in seq(along=model$comp) ) {
  
    if ( model$comp[[1]]$ngParam ) {
      scalePred <- sqrt(var(exp(model$comp[[j]]$f)))
    } else {
      scalePred <- sqrt(var(model$comp[[j]]$f)) 
    }
  
    ## Info from Martino Paper:
    ## 'True f' from Figure 3.
    truef <- c(0, 1.6, 2.6, 2.5, 2.6, 1.6, 0.9)
    truef <- truef/sqrt(var(truef))*scalePred
    
    ## Figure 2(a) histograms;
    B <- c(2.6, 1.5, 0.5, 0.2, 1.35) ## From Martino paper ... but don't know the scale
    B <- B/mean(B)*mean(modelB)      ## do a rough rescaling so that the scales match.
    S <- c(3, 0.8, 0.7, 1.8, 0.7)/1.8  
    D <- c(1.2, 1.6, 1.75, 3.2, 2.3)*0.8/3.2; 
  
    ## Martino f from Figure 2(b), again measured with a ruler.
    barencof <- t(matrix(c(0.0000000, 200.5201100, 355.5216125, 205.7574913, 135.0911372, 
              145.1080997, 130.7046969, 0.0000000, 184.0994134, 308.4759200, 
              232.1775328, 153.6595161, 85.7272235,  168.0910562, 0.0000000, 
              230.2262511, 337.5994811, 276.9416540, 164.5044287, 127.8653452, 173.6112139), 7, 3))
    barencof <- barencof/(1.8*mean(S))*mean(modelS)
    barencof <- barencof/(sqrt(apply(barencof, 1, var))%*%matrix(1,1,7))*scalePred

    if ( nonLinearity == "linear" ) {
      fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
      fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
      plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
    } else if ( nonLinearity == "multi" ) {
      fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
      fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
      plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
    } else {
      funcName <- model$comp[[j]]$nonLinearity
      func <- get(funcName, mode="function")
      fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
      fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
      plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
      title("Predicted Protein Concentration for p53")
      lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
    }

    for ( index in seq(length=model$comp[[j]]$numGenes) ) {
      xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]))+0.1
      xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]))-0.1
      plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
      title(paste(index, "-th Gene", sep=""))
      lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
      lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
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



gpsimMapResults <- function (model, expType, expNo, genes, scale=1, option=1) {
  if ( model$comp[[1]]$ngParam ) {
    nonLinearity <- "multi"
  } else {
    nonLinearity <- model$comp[[1]]$nonLinearity
  }

  if ( option == 1 ) {

    for ( j in seq(along=model$comp) ) {
      x11()
      if ( nonLinearity == "linear" ) {
        fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
        fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
        plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      } else if ( nonLinearity == "multi" ) {
        fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      } else {
        funcName <- model$comp[[j]]$nonLinearity
        func <- get(funcName, mode="function")
        fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      }

      uci <- model$comp[[j]]$y+2*sqrt(model$comp[[j]]$yvar)
      lci <- model$comp[[j]]$y-2*sqrt(model$comp[[j]]$yvar)                                      
      for ( index in seq(length=model$comp[[j]]$numGenes) ) {
        x11()
        xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), uci[,index])+0.1
        xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lci[,index])-0.1
        plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
        title(paste("Replica ", j, ": ", index, "-th Gene", genes[index], sep=""))
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        points(model$comp[[j]]$t, model$comp[[j]]$y[,index], lwd=3, col=3)
        arrows(model$comp[[j]]$t, lci[, index], model$comp[[j]]$t, uci[, index], length = 0.1, angle = 90, code = 3, col=3)
      }    
      
    }

    par(mar = c(7, 5, 6, 4))
    x11()
    bp <- barplot(model$comp[[1]]$B, mfg=c(3,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Basal Transcription Rates")
    x11()
    bp <- barplot(model$comp[[1]]$D, mfg=c(2,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Decay rate of the mRNA")
    x11()
    bp <- barplot(model$comp[[1]]$S, mfg=c(1,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Sensitivities to the Transcription Factor")
  } else {
 
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
    
    for ( j in seq(along=model$comp) ) {
      if ( nonLinearity == "linear" ) {
        fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
        fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
        plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      } else if ( nonLinearity == "multi" ) {
        fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      } else {
        funcName <- model$comp[[j]]$nonLinearity
        func <- get(funcName, mode="function")
        fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration")
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      }

      uci <- model$comp[[j]]$y+2*sqrt(model$comp[[j]]$yvar)
      lci <- model$comp[[j]]$y-2*sqrt(model$comp[[j]]$yvar) 
      for ( index in seq(length=model$comp[[j]]$numGenes) ) {
        xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), uci[,index])+0.1
        xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lci[,index])-0.1
        plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
        title(paste(index, "-th Gene", genes[index], sep=""))
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)

        points(model$comp[[j]]$t, model$comp[[j]]$y[,index], lwd=3, col=3)
        arrows(model$comp[[j]]$t, lci[, index], model$comp[[j]]$t, uci[, index], length = 0.1, angle = 90, code = 3, col=3)
      }    
      
    }

    par(mar = c(7, 5, 6, 4))
    bp <- barplot(model$comp[[1]]$B, mfg=c(3,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Basal Transcription Rates")
    bp <- barplot(model$comp[[1]]$D, mfg=c(2,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Decay rate of the mRNA")
    bp <- barplot(model$comp[[1]]$S, mfg=c(1,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Sensitivities to the Transcription Factor")
     
    dev.off()
    
  }
}



gpsimMapElkResults <- function (model, expType, expNo, genes, scale=1, option=1) {
  if ( model$comp[[1]]$ngParam ) {
    nonLinearity <- "multi"
  } else {
    nonLinearity <- model$comp[[1]]$nonLinearity
  }

  if ( option == 1 ) {

    for ( j in seq(along=model$comp) ) {
      x11()
      if ( nonLinearity == "linear" ) {
        fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
        fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
        plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      } else if ( nonLinearity == "multi" ) {
        fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      } else {
        funcName <- model$comp[[j]]$nonLinearity
        func <- get(funcName, mode="function")
        fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      }

      uci <- model$comp[[j]]$y+2*sqrt(model$comp[[j]]$yvar)
      lci <- model$comp[[j]]$y-2*sqrt(model$comp[[j]]$yvar)                                      
      for ( index in seq(length=model$comp[[j]]$numGenes) ) {
        x11()
        xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), uci[,index])+0.1
        xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lci[,index])-0.1
        plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
        title(paste("Replica ", j, ": ", index, "-th Gene", genes[index], sep=""))
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        points(model$comp[[j]]$t, model$comp[[j]]$y[,index], lwd=3, col=3)
        arrows(model$comp[[j]]$t, lci[, index], model$comp[[j]]$t, uci[, index], length = 0.1, angle = 90, code = 3, col=3)
      }    
      
    }

    par(mar = c(7, 5, 6, 4))
    x11()
    bp <- barplot(model$comp[[1]]$B, mfg=c(3,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Basal Transcription Rates")
    x11()
    bp <- barplot(model$comp[[1]]$D, mfg=c(2,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Decay rate of the mRNA")
    x11()
    bp <- barplot(model$comp[[1]]$S, mfg=c(1,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Sensitivities to the Transcription Factor")
  } else {
 
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
    
    for ( j in seq(along=model$comp) ) {
      if ( nonLinearity == "linear" ) {
        fmax <- max(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf))
        fmin <- min(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf))
        plot(model$comp[[j]]$mapt, model$comp[[j]]$f, type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf), lty=2, lwd=3, col=2)
      } else if ( nonLinearity == "multi" ) {
        fmax <- max(exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, exp(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, exp(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      } else {
        funcName <- model$comp[[j]]$nonLinearity
        func <- get(funcName, mode="function")
        fmax <- max(func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)))
        fmin <- min(func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)))
        plot(model$comp[[j]]$mapt, func(model$comp[[j]]$f), type="l", ylim=c(fmin,fmax), lwd=3, xlab="Time",ylab="")
        title("Predicted Protein Concentration for p53")
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f+2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, func(model$comp[[j]]$f-2*sqrt(model$comp[[j]]$varf)), lty=2, lwd=3, col=2)
      }

      uci <- model$comp[[j]]$y+2*sqrt(model$comp[[j]]$yvar)
      lci <- model$comp[[j]]$y-2*sqrt(model$comp[[j]]$yvar) 
      for ( index in seq(length=model$comp[[j]]$numGenes) ) {
        xmax <- max(model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), uci[,index])+0.1
        xmin <- min(model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lci[,index])-0.1
        plot(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index], type="l", ylim=c(xmin,xmax), lwd=3, xlab="Time",ylab="")
        title(paste(index, "-th Gene", genes[index], sep=""))
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]+2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)
        lines(model$comp[[j]]$mapt, model$comp[[j]]$ypred[,index]-2*sqrt(model$comp[[j]]$ypredVar[,index]), lty=2, lwd=3, col=2)

        points(model$comp[[j]]$t, model$comp[[j]]$y[,index], lwd=3, col=3)
        arrows(model$comp[[j]]$t, lci[, index], model$comp[[j]]$t, uci[, index], length = 0.1, angle = 90, code = 3, col=3)
      }    
      
    }

    par(mar = c(7, 5, 6, 4))
    bp <- barplot(model$comp[[1]]$B, mfg=c(3,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Basal Transcription Rates")
    bp <- barplot(model$comp[[1]]$D, mfg=c(2,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Decay rate of the mRNA")
    bp <- barplot(model$comp[[1]]$S, mfg=c(1,1))
    text(bp, par("usr")[3]-0.05, srt=45, adj=1, labels=genes, xpd=TRUE)
    title("Sensitivities to the Transcription Factor")
     
    dev.off()
    
  }
}



gpsimMapInitParam <- function (y, yvar, options, display=TRUE, sopt=FALSE) {
  numGenes <- dim(y)[2]
  if ( sopt ) {
    optionSset <- c(0.01, 1, 10, 100, 1000)
  } else {
    if ( (!any(names(options)=="S")) ) {
      optionSset <- 1
    } else {
      optionSset <- options$S[1]
    }
  }
  
  initOptions <- options     
  mu <- apply(y, 2, mean)
  optionDset <- c(0.01, 0.1, 0.4, 1)
  set.seed(10)
  model <- list(gene=list())
  
  for ( j in seq(length=numGenes) ) {
    initOptions$gParam <- options$gParam[j]     
    model$gene[[j]] <- list()
    if ( display ) {
      x11()
      ymin = -1
      ymax = max(y[,j])*2
      plot(options$times, y[,j], ylim=c(ymin, ymax))
    }

    for ( si in seq(along=optionSset) ) {
      initOptions$S <- optionSset[si]
      for ( dInd in seq(along=optionDset) ) {
        di <- (si-1)*length(optionSset) + dInd
        initOptions$D <- optionDset[dInd]
        initOptions$B <- initOptions$D*mu[j]
        
        model$gene[[j]]$init[[di]] <- gpsimMapCreate(1, 1, options$times, as.matrix(y[,j]), as.matrix(yvar[,j]), initOptions)
        if ( length(options$kern)>1 ) {
          if ( options$kern[2]=="mlp" ) {
            model$gene[[j]]$init[[di]]$kern$weightVariance <- 30
            model$gene[[j]]$init[[di]]$kern$biasVariance <- 20
            params <- gpsimMapExtractParam(model$gene[[j]]$init[[di]])
            model$gene[[j]]$init[[di]] <- gpsimMapExpandParam(model$gene[[j]]$init[[di]], params)
          } else if ( options$kern[2]=="rbf" ) {
            model$gene[[j]]$init[[di]]$kern$inverseWidth <- 0.001
            params <- gpsimMapExtractParam(model$gene[[j]]$init[[di]])
            model$gene[[j]]$init[[di]] <- gpsimMapExpandParam(model$gene[[j]]$init[[di]], params)
          }
        } else {
          if ( options$kern=="mlp" ) {
            model$gene[[j]]$init[[di]]$kern$weightVariance <- 30
            model$gene[[j]]$init[[di]]$kern$biasVariance <- 20
            params <- gpsimMapExtractParam(model$gene[[j]]$init[[di]])
            model$gene[[j]]$init[[di]] <- gpsimMapExpandParam(model$gene[[j]]$init[[di]], params)
          } else if ( options$kern=="rbf" ) {
            model$gene[[j]]$init[[di]]$kern$inverseWidth <- 0.001
            params <- gpsimMapExtractParam(model$gene[[j]]$init[[di]])
            model$gene[[j]]$init[[di]] <- gpsimMapExpandParam(model$gene[[j]]$init[[di]], params)
          }
        }
        
        updateFoptions <- optimiFdefaultOptions()
        model$gene[[j]]$init[[di]] <- gpsimMapUpdateF(model$gene[[j]]$init[[di]], updateFoptions)
        
        ypredMean <- mean(model$gene[[j]]$init[[di]]$ypred)
        ypredScale <- sqrt(var(model$gene[[j]]$init[[di]]$ypred[model$gene[[j]]$init[[di]]$timesIndex]))
        model$gene[[j]]$init[[di]]$B <- (model$gene[[j]]$init[[di]]$B - initOptions$D*(ypredMean - ypredScale*mean(model$gene[[j]]$init[[di]]$y)))/ypredScale

        if ( model$gene[[j]]$init[[di]]$B < 0 && any(grep("isGroupNonlinearity", names(model$gene[[j]]$init[[di]]))) ) {
       #   if ( model$gene[[j]]$init[[di]]$includeRepression ) {
            model$gene[[j]]$init[[di]]$alpha <- model$gene[[j]]$init[[di]]$alpha/ypredScale
            model$gene[[j]]$init[[di]]$alpha <- model$gene[[j]]$init[[di]]$alpha + model$gene[[j]]$init[[di]]$B
            model$gene[[j]]$init[[di]]$B <- 1e-6
       #   } else {
       #     model$gene[[j]]$init[[di]]$B <- 1e-6
       #   }
        }
        
        model$gene[[j]]$init[[di]]$S <- model$gene[[j]]$init[[di]]$S/ypredScale 
        params <- gpsimMapExtractParam(model$gene[[j]]$init[[di]])  
        model$gene[[j]]$init[[di]] <- gpsimMapExpandParam(model$gene[[j]]$init[[di]], params)
        f <- gpsimMapFunctionalExtractParam(model$gene[[j]]$init[[di]])
        model$gene[[j]]$init[[di]] <- gpsimMapFunctionalExpandParam(model$gene[[j]]$init[[di]], f)
        model$gene[[j]]$ll[di] <- gpsimMapLogLikelihood(model$gene[[j]]$init[[di]])
        if ( display )
          lines(model$gene[[j]]$init[[di]]$mapt, model$gene[[j]]$init[[di]]$ypred)
      }
    }
   
    index <- order(model$gene[[j]]$ll, decreasing=TRUE)[1]
    if ( display )
      lines(model$gene[[j]]$init[[index]]$mapt, model$gene[[j]]$init[[index]]$ypred, col=2)    
    
    options$D[j] <- model$gene[[j]]$init[[index]]$D
    options$S[j] <- model$gene[[j]]$init[[index]]$S

    if ( any(names(model$gene[[j]]$init[[index]])=="alpha") ) {
      if ( j == 1 )
        options$alpha <- array()
      options$alpha[j] <- model$gene[[j]]$init[[index]]$alpha
    }
    options$B[j] <- model$gene[[j]]$init[[index]]$B     
    options$gParam[j] <- model$gene[[j]]$init[[index]]$gParam
  }

  set.seed(11)
  if ( !sopt ) 
    options$S <- runif(numGenes)
  options$B <- options$D*mu

  return(options)

}
