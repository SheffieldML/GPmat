simKernParamInit <- function (kern) {

  if ( kern$inputDimension > 1 )
    stop("SIM kernel is only valid for one-D input.")
  kern$delay <- 0
  kern$decay <- 1
  kern$initVal <- 1
  kern$variance <- 1
  kern$inverseWidth <- 1
  if ("options" %in% names(kern) && "gaussianInitial" %in% names(kern$options) && kern$options$gaussianInitial) {
    kern$gaussianInitial <- TRUE
    kern$initialVariance <- 1
    kern$nParams <- 4
    kern$paramNames <- c("decay", "inverseWidth", "variance", "initialVariance")
  }
  else {
    kern$gaussianInitial <- FALSE
    kern$nParams <- 3
    kern$paramNames <- c("decay", "inverseWidth", "variance")
  }

  if ("options" %in% names(kern) && "isNegativeS" %in% names(kern$options) && kern$options$isNegativeS) {
    kern$isNegativeS <- TRUE
    kern$transforms <- list(index=setdiff(1:kern$nParams, 3), type="positive")
    kern$paramNames[3] <- "sensitivity"
  }
  else {
    kern$isNegativeS <- FALSE
    kern$transforms <- list(index=1:kern$nParams, type="positive")
  }
  kern$isStationary <- FALSE

  return (kern)

}

simXrbfKernCompute <- function (simKern, rbfKern, t1, t2=t1) {
  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if (rbfKern$variance != 1)
    warning('Non-unit RBF kernel variance. The SIM kernel can only be cross combined with an RBF kernel with variance 1.')

  if ( simKern$inverseWidth != rbfKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1 <- t1 - simKern$delay
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat
  sigma <- sqrt(2/simKern$inverseWidth)

  invSigmaDiffT <- 1/sigma*diffT
  halfSigmaDi <- 0.5*sigma*simKern$decay

  lnPart <- lnDiffErfs(halfSigmaDi + t2Mat/sigma, halfSigmaDi - invSigmaDiffT)

  K <- lnPart[[2]] * exp(halfSigmaDi*halfSigmaDi - simKern$decay*diffT + lnPart[[1]])

  #K <- 0.5*sqrt(simKern$variance)*sqrt(rbfKern$variance)*K*sqrt(pi)*sigma
  K <- K*sqrt(pi)
#  if(kern$isNegativeS)
    K <- 0.5*simKern$sensitivity*K*sigma
#  else
#    K <- 0.5*sqrt(simKern$variance)*K*sigma
  return (K)
}



simXsimKernCompute <- function (simKern1, simKern2, t1, t2=t1) {
  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if ( simKern1$inverseWidth != simKern2$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  sigma <- sqrt(2/simKern1$inverseWidth)
  h1 <- simComputeH(t1, t2, simKern1$decay, simKern2$decay, simKern1$delay, simKern2$delay, sigma)
  h2 <- simComputeH(t2, t1, simKern2$decay, simKern1$decay, simKern2$delay, simKern1$delay, sigma)
  K <- h1 + t(h2)
  K <- 0.5*K*sqrt(pi)*sigma
#  if(simKern1$isNegativeS)
  K <- simKern1$sensitivity*simKern2$sensitivity*K
#  else
#    K <- sqrt(simKern1$variance)*K
#  if(simKern2$isNegativeS)
#  K <- simKern2$sensitivity*K
#  else
#    K <- sqrt(simKern2$variance)*K

  return (K)
}



simKernExtractParam <- function (kern, option=1) {

  if (kern$gaussianInitial) {
    if (kern$isNegativeS) {
      if ( option == 1 ) {
        params <- c(kern$decay, kern$inverseWidth, kern$sensitivity, kern$initialVariance)
      } else {
        params <- list(values=c(kern$decay, kern$inverseWidth, kern$sensitivity, kern$initialVariance), names=kern$paramNames)
      }
    }
    else {
      if ( option == 1 ) {
        params <- c(kern$decay, kern$inverseWidth, kern$variance, kern$initialVariance)
      } else {
        params <- list(values=c(kern$decay, kern$inverseWidth, kern$variance, kern$initialVariance), names=kern$paramNames)
      }
    }
  }
  else {
    if (kern$isNegativeS) {
      if ( option == 1 ) {
        params <- c(kern$decay, kern$inverseWidth, kern$sensitivity)
      } else {
        params <- list(values=c(kern$decay, kern$inverseWidth, kern$sensitivity), names=kern$paramNames)
      }
    }
    else {
      if ( option == 1 ) {
        params <- c(kern$decay, kern$inverseWidth, kern$variance)
      } else {
        params <- list(values=c(kern$decay, kern$inverseWidth, kern$variance), names=kern$paramNames)
      }
    }
  }
  return (params)
}



simKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$decay <- params[1]
  kern$inverseWidth <- params[2]
  if(kern$isNegativeS) {
    kern$sensitivity <- params[3]
    kern$variance <- kern$sensitivity*kern$sensitivity
  }
  else {
    kern$variance <- params[3]
    kern$sensitivity <- sqrt(kern$variance)
  }
  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial )
    kern$initialVariance <- params[4]

  return (kern)
}

simKernDisplay <- function (kern, spaceNum=0) {

  spacing = matrix("", spaceNum+1)
  cat(spacing)
  if(kern$isStationary)
    cat("Stationary version of the kernel\n")
  else
    cat("Non-stationary version of the kernel\n")

##   cat(spacing)
##   if(kern$isNormalised)
##     cat("Normalised version of the kernel\n")
##   else
##     cat("Unnormalised version of the kernel\n")

  if(kern$isNegativeS)
    cat("Sensitivities allowed to be negative.\n")
  else
    cat("Sensitivities constrained positive.\n")
  cat(spacing)
  cat("SIM decay: ", kern$decay, "\n", sep="")
  cat(spacing)
  cat("SIM inverse width: ", kern$inverseWidth, " (length scale ", 1/sqrt(kern$inverseWidth), ")\n", sep="")
  cat(spacing)
  if(kern$isNegativeS) 
    cat("SIM Sensitivity: ", kern$sensitivity, "\n", sep="")  
  else
    cat("SIM Variance: ", kern$variance, "\n", sep="")
  if(kern$gaussianInitial) {
    cat(spacing)
    cat("SIM Initial Variance: ", kern$initialVariance, "\n", sep="")
  }
  ## cat(spacing)
  ## cat("SIM delay: %2.4f\n", kern$delay)
}


simKernCompute <- function (kern, t, t2=t) {
  if ( ( dim(as.matrix(t))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  sigma <- sqrt(2/kern$inverseWidth)

  h <- simComputeH(t, t2, kern$decay, kern$decay, kern$delay, kern$delay, sigma)

  if ( nargs() < 3 ) {
    k <- h + t(h)
  } else {
    h2 <- simComputeH(t2, t, kern$decay, kern$decay, kern$delay, kern$delay, sigma)
    k <- h + t(h2)
  }

  k <- 0.5*k*sqrt(pi)*sigma
  k <- kern$variance*k

  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial ) {
    dim1 <- dim(as.matrix(t))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    k = k + kern$initialVariance * exp(- kern$decay * (t1Mat + t2Mat))
  }

  return (k)
}



simComputeH <-  function (t1, t2, Di, Dj, deltai, deltaj, sigma, option=1) {
  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1 <- t1 - deltai
  t2 <- t2 - deltaj
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat
  invSigmaDiffT <- 1/sigma*diffT
  halfSigmaDi <- 0.5*sigma*Di
  h <- matrix(0, dim1, dim2)

  lnPart1_full <- lnDiffErfs(halfSigmaDi + t2Mat/sigma, halfSigmaDi - invSigmaDiffT)
  lnPart2_full <- lnDiffErfs(halfSigmaDi, halfSigmaDi - t1Mat/sigma)

  lnPart1 <- lnPart1_full[[1]]
  signs1 <- lnPart1_full[[2]]
  lnPart2 <- lnPart2_full[[1]]
  signs2 <- lnPart2_full[[2]]

  h <- signs1 * exp(halfSigmaDi*halfSigmaDi-Di*diffT+lnPart1-log(Di+Dj)) - signs2 * exp(halfSigmaDi*halfSigmaDi-Di*t1Mat-Dj*t2Mat+lnPart2-log(Di + Dj))

  sigma2 <- sigma*sigma

  if ( option == 1 ) {
    return (h)
  } else {
    dhdDi <- (0.5*Di*sigma2*(Di+Dj)-1)*h + (-diffT*signs1*exp(halfSigmaDi*halfSigmaDi-Di*diffT+lnPart1) + t1Mat*signs2*exp(halfSigmaDi*halfSigmaDi-Di*t1Mat-Dj*t2Mat+lnPart2)) + sigma/sqrt(pi)*(-exp(-diffT*diffT/sigma2)+exp(-t2Mat*t2Mat/sigma2-Di*t1Mat)+exp(-t1Mat*t1Mat/sigma2-Dj*t2Mat)-exp(-(Di*t1Mat+Dj*t2Mat)))
    dhdDi <- Re(dhdDi/(Di+Dj))

    dhdDj <- t2Mat*signs2*exp(halfSigmaDi*halfSigmaDi-(Di*t1Mat+Dj*t2Mat)+lnPart2)-h
    dhdDj <- Re(dhdDj/(Di+Dj))

    dhdsigma <- 0.5*Di*Di*sigma*h + 2/(sqrt(pi)*(Di+Dj))*((-diffT/sigma2-Di/2)*exp(-diffT*diffT/sigma2) + (-t2Mat/sigma2+Di/2)*exp(-t2Mat*t2Mat/sigma2-Di*t1Mat) - (-t1Mat/sigma2-Di/2)*exp(-t1Mat*t1Mat/sigma2-Dj*t2Mat) - Di/2*exp(-(Di*t1Mat+Dj*t2Mat)))

    simH <- list(h=h, dhdDi=dhdDi, dhdDj=dhdDj, dhdsigma=dhdsigma)
    return (simH)
  }
}



simKernGradient <- function (kern, t, t2, covGrad) {
  if ( nargs() == 3 ) {
    covGrad <- t2
    t2 <- t
  }

  gFull <- simXsimKernGradient(kern, kern, t, t2, covGrad)

  g <- gFull$g1 + gFull$g2

  return (g)
}



simXsimKernGradient <- function (simKern1, simKern2, t1, t2, covGrad) {
  if ( nargs() < 5 ) {
    covGrad <- t2
    t2 <- t1
  }

  if ( dim(as.matrix(t1))[2]>1 | dim(as.matrix(t2))[2]>1 )
    stop("Input can only have one column for SIM kernel.")

  if ( simKern1$inverseWidth != simKern2$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  option <- 2
  sigma <- sqrt(2/simKern1$inverseWidth)
  simH1 <- simComputeH(t1, t2, simKern1$decay, simKern2$decay, simKern1$delay, simKern2$delay, sigma, option)
  h1 <- simH1$h
  dh1dD1 <- simH1$dhdDi
  dh1dD2 <- simH1$dhdDj
  dh1dsigma <- simH1$dhdsigma

  simH2 <- simComputeH(t2, t1, simKern2$decay, simKern1$decay, simKern2$delay, simKern1$delay, sigma, option)
  h2 <- simH2$h
  dh2dD2 <- simH2$dhdDi
  dh2dD1 <- simH2$dhdDj
  dh2dsigma <- simH2$dhdsigma

  dKdD1 <- dh1dD1 + t(dh2dD1)
  dKdD2 <- dh1dD2 + t(dh2dD2)
  dKdsigma <- dh1dsigma + t(dh2dsigma)

  C1 <- simKern1$sensitivity
  C2 <- simKern2$sensitivity

  K <- h1 + t(h2)
  K <- 0.5*K*sqrt(pi)

  var2 <- C1*C2
  dkdD1 <- (sum(covGrad*dh1dD1)+sum(covGrad*t(dh2dD1)))*0.5*sqrt(pi)*sigma*var2
  dkdD2 <- (sum(covGrad*dh1dD2)+sum(covGrad*t(dh2dD2)))*0.5*sqrt(pi)*sigma*var2
  dkdsigma <- sum(covGrad*(dKdsigma*0.5*sqrt(pi)*sigma + K))*var2
  K <- sigma*K
  dkdC1 <- C2*sum(covGrad*K)
  dkdC2 <- C1*sum(covGrad*K)

  if(simKern1$isNegativeS)
    dkdSim1Variance <- dkdC1
  else
    dkdSim1Variance <- dkdC1*0.5/C1

  if(simKern2$isNegativeS)
    dkdSim2Variance <- dkdC2
  else
    dkdSim2Variance <- dkdC2*0.5/C2

  dkdinvWidth <- -0.5*sqrt(2)/(simKern1$inverseWidth*sqrt(simKern1$inverseWidth))*dkdsigma
  K <- var2*K

  g1 <- c(dkdD1, dkdinvWidth, dkdSim1Variance)
  g2 <- c(dkdD2, 0, dkdSim2Variance)

  g <- list(g1=g1, g2=g2)
  return (g)
}



simXrbfKernGradient <- function (simKern, rbfKern, t1, t2, covGrad) {
  if ( nargs() < 5 ) {
    covGrad <- t2
    t2 <- t1
  }

  if ( dim(as.matrix(t1))[2]>1 | dim(as.matrix(t2))[2]>1 )
    stop("Input can only have one column for SIM kernel.")

  if ( simKern$inverseWidth != rbfKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( nargs()<5 ) {
    k <- simXrbfKernCompute(simKern, rbfKern, t1)
  } else {
    k <- simXrbfKernCompute(simKern, rbfKern, t1, t2)
  }

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat
  sigma <- sqrt(2/simKern$inverseWidth)
  sigma2 <- sigma*sigma
  Ci <- simKern$sensitivity
  Cj <- sqrt(rbfKern$variance)
  Di <- simKern$decay

  part2 <- exp(-t2Mat*t2Mat/sigma2-t1Mat*Di) - exp(-diffT*diffT/sigma2)

  dkdD <- sum((k*(0.5*sigma2*Di - diffT) + 0.5*Ci*Cj*sigma2*part2)*covGrad)

  dkdsigma <- sum((k*(1/sigma+0.5*sigma*Di*Di)+Ci*Cj*sigma*((-diffT/sigma2-Di/2)*exp(-diffT*diffT/sigma2)+(-t2Mat/sigma2+Di/2)*exp(-t2Mat*t2Mat/sigma2-t1Mat*Di)))*covGrad)

  dkdC <- sum(k*covGrad)/Ci
  dkdRbfVariance <- 0.5*sum(k*covGrad)/rbfKern$variance

  dkdinvWidth <-  -0.5*sqrt(2)/(simKern$inverseWidth*sqrt(simKern$inverseWidth))*dkdsigma
  dkdSimVariance <- dkdC*0.5/Ci

  g1 <- c(dkdD, dkdinvWidth, dkdSimVariance)
  g2 <- c(0, dkdRbfVariance)

  g <- list(g1=g1, g2=g2)
  return (g)
}



simKernDiagCompute <- function (kern, t) {
  if ( dim(as.matrix(t))[2]>1 )
    stop("Input can only have one column.")

  sigma <- sqrt(2/kern$inverseWidth)
  t <- t - kern$delay
  halfSigmaD <- 0.5*sigma*kern$decay

  lnPart1_full <- lnDiffErfs(halfSigmaD + t/sigma, halfSigmaD)
  lnPart2_full <- lnDiffErfs(halfSigmaD, halfSigmaD - t/sigma)

  lnPart1 <- lnPart1_full[[1]]
  signs1 <- lnPart1_full[[2]]
  lnPart2 <- lnPart2_full[[1]]
  signs2 <- lnPart2_full[[2]]

  h <- signs1 * exp(halfSigmaD*halfSigmaD + lnPart1) - signs2 * exp(halfSigmaD*halfSigmaD-(2*kern$decay*t) + lnPart2)

  # h <- exp(halfSigmaD*halfSigmaD)*(erff(-halfSigmaD)+erff(t/sigma+halfSigmaD))-exp(halfSigmaD*halfSigmaD-(2*kern$decay*t))*(erff(t/sigma-halfSigmaD)+erff(halfSigmaD))

  k <- 2*h
  k <- 0.5*k*sqrt(pi)*sigma
  k <- kern$variance*k/(2*kern$decay)
  return (k)
}


