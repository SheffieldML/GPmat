disimKernParamInit <- function (kern) {

  if ( kern$inputDimension > 1 )
    stop("SIM kernel is only valid for one-D input.")

  kern$di_decay <- 0.1
  kern$inverseWidth <- 1
  kern$di_variance <- 1
  kern$decay <- 1
  kern$variance <- 1
  kern$rbf_variance <- 1

  if ("options" %in% names(kern) && "gaussianInitial" %in% names(kern$options) && kern$options$gaussianInitial) {
    kern$gaussianInitial <- TRUE
    kern$initialVariance <- 1
    kern$nParams <- 7
    kern$paramNames <- c("di_decay", "inverseWidth", "di_variance", "decay", "variance", "rbf_variance", "initialVariance")
  }
  else {
    kern$gaussianInitial <- FALSE
    kern$nParams <- 6
    kern$paramNames <- c("di_decay", "inverseWidth", "di_variance", "decay", "variance", "rbf_variance")
  }

  if ("options" %in% names(kern) && "isNormalised" %in% names(kern$options) && kern$options$isNormalised) {
    kern$isNormalised <- TRUE
  }
  else
    kern$isNormalised <- FALSE

  if ("options" %in% names(kern) && "inverseWidthBounds" %in% names(kern$options)) {
    kern$transforms <- list(list(index=setdiff(1:kern$nParams, 2), type="positive"),
                            list(index=2, type="bounded"))
    kern$transformArgs <- list()
    kern$transformArgs[[2]] <- kern$options$inverseWidthBounds
    kern$inverseWidth <- mean(kern$options$inverseWidthBounds)
  }
  else
    kern$transforms <- list(list(index=1:kern$nParams, type="positive"))

  kern$isStationary <- FALSE

  return (kern)
}



disimKernCompute <- function (kern, t1, t2=t1) {

  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  l <- sqrt(2/kern$inverseWidth)

  h <- .disimComputeH(t1, t2, kern$di_decay, kern$decay, kern$decay, l)
  hp <- .disimComputeHPrime(t1, t2, kern$di_decay, kern$decay, kern$decay, l)
  if ( nargs()<3 ) {
    k <- h + t(h) + hp + t(hp)
  } else {
    h2 <- .disimComputeH(t2, t1, kern$di_decay, kern$decay, kern$decay, l)
    hp2 <- .disimComputeHPrime(t2, t1, kern$di_decay, kern$decay, kern$decay, l)
    k <- h + t(h2) + hp + t(hp2)
  }
  k <- 0.5*kern$rbf_variance*kern$di_variance*kern$variance*k
  if (!kern$isNormalised)
    k <- k*sqrt(pi)*l
  k <- Re(k)

  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial ) {
    dim1 <- dim(as.matrix(t1))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t1, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    k = k + kern$initialVariance * kern$variance * 
      (exp(- kern$di_decay * t1Mat) - exp(- kern$decay * t1Mat)) /
        (kern$decay - kern$di_decay) * 
          (exp(- kern$di_decay * t2Mat) - exp(- kern$decay * t2Mat)) /
            (kern$decay - kern$di_decay);
  }

  return (k)
}



.disimComputeH <-  function (t1, t2, delta, Dj, Dk, l, option=1) {

  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  invLDiffT <- 1/l*diffT
  halfLDelta <- 0.5*l*delta
  h <- matrix(0, dim1, dim2)

  lnPart1_f <- lnDiffErfs(halfLDelta - t1Mat/l, halfLDelta)
  lnPart2_f <- lnDiffErfs(halfLDelta + t2Mat/l, halfLDelta-invLDiffT)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]
  lnPart2 <- lnPart2_f[[1]]
  signs2 <- lnPart2_f[[2]]

  lnCommon <- halfLDelta^2-log(2*delta)-Dk*t2Mat-delta*t1Mat-log(0i + Dj-delta)

  lnFact1a1 <- log(Dk + delta) + (Dk-delta)*t2Mat
  lnFact1a2 <- log(2*delta)
  lnFact1b <- log(0i + Dk^2 - delta^2)

  lnFact2a <- (Dk + delta)*t2Mat
  lnFact2b <- log(Dk + delta)

  if (abs(Dk - delta) < .1)
    h <- signs1 * exp(lnCommon + lnPart1) * ((exp((Dk-delta)*t2Mat) - 1) / (Dk - delta) + 1/(Dk+delta)) + signs2 * exp(lnCommon + lnFact2a - lnFact2b + lnPart2)
  else
    h <- signs1 * exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) - signs1 * exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) + signs2 * exp(lnCommon + lnFact2a - lnFact2b + lnPart2)

  h <- Re(h)

  l2 <- l*l

  if ( option == 1 ) {
    return (h)
  } else {
    m1 <- pmin((halfLDelta - t1Mat/l)^2, halfLDelta^2)
    m2 <- pmin((halfLDelta + t2Mat/l)^2, (halfLDelta - invLDiffT)^2)

    dlnPart1 <- l/sqrt(pi)*(exp(-(halfLDelta - t1Mat/l)^2 + m1) - exp(-halfLDelta^2 + m1))
    dlnPart2 <- l/sqrt(pi)*(exp(-(halfLDelta + t2Mat/l)^2 + m2) - exp(-(halfLDelta - invLDiffT)^2 + m2))
    dh_ddelta <- (l*halfLDelta - t1Mat + 1/(Dj-delta)-1/delta)*h + exp(lnCommon + lnFact1a1 - lnFact1b - m1)*dlnPart1 - exp(lnCommon + lnFact1a2 - lnFact1b - m1)*dlnPart1 + exp(lnCommon + lnFact2a - lnFact2b - m2)*dlnPart2 + signs1 * exp(lnCommon + lnPart1 + lnFact1a1 - lnFact1b)*(1/(Dk-delta) - t2Mat) - signs1 * exp(lnCommon + lnPart1 + log(2*(Dk^2 + delta^2)) -2*lnFact1b) + (t2Mat - 1/(Dk + delta))*signs2*exp((Dk + delta)*t2Mat + lnCommon + lnPart2 - log(Dk + delta))
    dh_ddelta <- Re(dh_ddelta)

    disimH <- list(h=h, dh_ddelta=dh_ddelta)

    if ( option > 2 ) {
      dh_dD_j = Re(- 1/(Dj - delta)*h)

      disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j)

      if ( option > 3 ) {
        dh_dD_k <- -t2Mat*h + signs1 * exp(lnCommon + lnPart1 + lnFact1a1 - lnFact1b) * (t2Mat - 1/(Dk - delta)) - signs1 * exp(lnCommon + lnPart1) * (-4*delta*Dk / (Dk^2 - delta^2)^2) + (t2Mat - 1/(Dk + delta)) * signs2 * exp(lnFact2a - lnFact2b + lnCommon + lnPart2)
        dh_dD_k <- Re(dh_dD_k)

        disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j, dh_dD_k=dh_dD_k)

        if ( option > 4 ) {
          dh_dl <- exp(lnCommon + lnFact1a1 - lnFact1b - m1) * ((delta/sqrt(pi) + 2*t1Mat/(l2*sqrt(pi))) * exp(-(halfLDelta - t1Mat/l)^2 + m1) - (delta/sqrt(pi) * exp(-halfLDelta^2 + m1))) - exp(lnCommon + lnFact1a2 - lnFact1b - m1) * ((delta/sqrt(pi) + 2*t1Mat/(l2*sqrt(pi))) * exp(-(halfLDelta - t1Mat/l)^2 + m1) - (delta/sqrt(pi) * exp(-halfLDelta^2 + m1))) +exp(lnCommon + lnFact2a - lnFact2b - m2) * ((delta/sqrt(pi) - 2*t2Mat/(l2*sqrt(pi))) * exp(-(halfLDelta + t2Mat/l)^2 + m2) - ((delta/sqrt(pi) + 2*invLDiffT/(l*sqrt(pi))) * exp(-(halfLDelta - invLDiffT)^2 + m2))) + delta*halfLDelta*h
          dh_dl <- Re(dh_dl)
          disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j, dh_dD_k=dh_dD_k, dh_dl=dh_dl)
        }
      }
    }
    return (disimH)
  }
}



.disimComputeHPrime <-  function (t1, t2, delta, Dj, Dk, l, option=1) {

  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  invLDiffT <- 1/l*diffT
  halfLD_k <- 0.5*l*Dk
  h <- matrix(0, dim1, dim2)

  lnPart1_f <- lnDiffErfs(halfLD_k - t2Mat/l, halfLD_k)
  lnPart2_f <- lnDiffErfs(halfLD_k + t1Mat/l, halfLD_k+invLDiffT)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]
  lnPart2 <- lnPart2_f[[1]]
  signs2 <- lnPart2_f[[2]]

  lnCommon <- halfLD_k^2 - Dj*t1Mat - Dk*t2Mat - log(0i + delta^2-Dk^2) - log(Dj+Dk)

  lnFact1a1 <- log(Dk + delta)
  lnFact1a2 <- log(Dj + Dk) + (Dj-delta)*t1Mat
  lnFact1b <- log(0i + delta-Dj)

  lnFact2a <- (Dj + Dk) * t1Mat

  if (abs(Dj - delta) < .1)
    h <- signs1 * exp(lnCommon + lnPart1) * (1 + (Dj+Dk)*(1-exp((Dj-delta)*t1Mat))/(delta-Dj)) + signs2 * exp(lnCommon + lnFact2a + lnPart2)
  else
    h <- signs1 * exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) - signs1 * exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) + signs2 * exp(lnCommon + lnFact2a + lnPart2)

  h <- Re(h)

  l2 <- l*l

  if ( option == 1 ) {
    return (h)
  } else {
    dh_ddelta <- -2*delta/(delta^2-Dk^2)*h - signs1 * exp(lnCommon + lnPart1 + log(Dk + Dj) - 2*lnFact1b) - signs1 * exp(lnCommon + lnPart1 + lnFact1a2 - lnFact1b)*(-t1Mat - 1/(delta - Dj))
    dh_ddelta <- Re(dh_ddelta)
    disimH <- list(h=h, dh_ddelta=dh_ddelta)

    if ( option > 2 ) {
      dh_dD_j <- (-t1Mat - 1 /(Dj + Dk))*h + signs1 * exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1)*(1/(delta - Dj)) - signs1 * exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1)*(1/(Dj + Dk) + t1Mat + 1/(delta - Dj)) + t1Mat*signs2*exp(lnCommon + lnPart2 + lnFact2a)
      dh_dD_j <- Re(dh_dD_j)

      disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j)
      if ( option > 3 ) {
        m1 <- pmin((halfLD_k - t2Mat/l)^2, halfLD_k^2)
        m2 <- pmin((halfLD_k + t1Mat/l)^2, (halfLD_k + invLDiffT)^2)

        dlnPart1 <- l/sqrt(pi) * (exp(-(halfLD_k - t2Mat/l)^2 + m1) - exp(-halfLD_k^2 + m1))
        dlnPart2 <- l/sqrt(pi) * (exp(-(halfLD_k + t1Mat/l)^2 + m2) - exp(-(halfLD_k + invLDiffT)^2 + m2))
        dh_dD_k <- (l*halfLD_k + 2*Dk / (delta^2 - Dk^2) + (-t2Mat - 1/(Dj + Dk)))*h + exp(lnCommon + lnFact1a1 - lnFact1b - m1) * dlnPart1 - exp(lnCommon + lnFact1a2 - lnFact1b - m1) * dlnPart1 + exp(lnCommon + lnFact2a - m2) * dlnPart2 + signs1 * exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) * (1/(Dk + delta)) - signs1 * exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) *(1/(Dj + Dk)) + signs2 * exp(lnCommon + lnPart2 + log(0i + t1Mat) + lnFact2a)
        dh_dD_k <- Re(dh_dD_k)

        disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j, dh_dD_k=dh_dD_k)

        if ( option > 4 ) {
          dh_dl <- exp(lnCommon + lnFact1a1 - lnFact1b - m1) * ((Dk/sqrt(pi) + 2*t2Mat/(l2*sqrt(pi))) * exp(-(halfLD_k - t2Mat/l)^2 + m1) - (Dk/sqrt(pi) * exp(-halfLD_k^2 + m1))) -exp(lnCommon + lnFact1a2 - lnFact1b - m1) * ((Dk/sqrt(pi) + 2*t2Mat/(l2*sqrt(pi))) * exp(-(halfLD_k - t2Mat/l)^2 + m1) - (Dk/sqrt(pi) * exp(-halfLD_k^2 + m1))) +exp(lnCommon + lnFact2a - m2) * ((Dk/sqrt(pi) - 2*t1Mat/(l2*sqrt(pi))) * exp(-(halfLD_k + t1Mat/l)^2 + m2) - ((Dk/sqrt(pi) - 2*invLDiffT/(l*sqrt(pi))) * exp(-(halfLD_k + invLDiffT)^2 + m2)))+ Dk*halfLD_k*h

          dh_dl <- Re(dh_dl)
          disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j, dh_dD_k=dh_dD_k, dh_dl=dh_dl)
        }
      }
    }
    return (disimH)
  }
}



disimKernExtractParam <- function (kern, only.values=TRUE) {
  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial )
    params <- c(kern$di_decay, kern$inverseWidth, kern$di_variance, kern$decay, kern$variance, kern$rbf_variance, kern$initialVariance)
  else
    params <- c(kern$di_decay, kern$inverseWidth, kern$di_variance, kern$decay, kern$variance, kern$rbf_variance)

  if ( !only.values )
    names(params) <- kern$paramNames
  
  return (params)
}



disimKernExpandParam <- function (kern, params) {
  if ( is.list(params) )
    params <- params$values

  kern$di_decay <- params[1]
  kern$inverseWidth <- params[2]
  kern$di_variance <- params[3]
  kern$decay <- params[4]
  kern$variance <- params[5]
  kern$rbf_variance <- params[6]
  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial )
    kern$initialVariance <- params[7]

  return (kern)
}


disimKernDisplay <- function (kern, spaceNum=0) {

  spacing = matrix("", spaceNum+1)
##   cat(spacing)
##   if(kern$isStationary)
##     cat("Stationary version of the kernel\n")
##   else
##     cat("Non-stationary version of the kernel\n")

  cat(spacing)
  if(kern$isNormalised)
    cat("Normalised version of the kernel\n")
  else
    cat("Unnormalised version of the kernel\n")

##   cat(spacing)
##   if(kern$isNegativeS)
##     cat("Sensitivities allowed to be negative.\n")
##   else
##     cat("Sensitivities constrained positive.\n")
  cat(spacing)
  cat("DISIM decay: ", kern$di_decay, "\n", sep="")
  cat(spacing)
  cat("DISIM inverse width: ", kern$inverseWidth, " (length scale ", 1/sqrt(kern$inverseWidth), ")\n", sep="")
  cat(spacing)
  cat("DISIM Variance: ", kern$di_variance, "\n", sep="")
  cat(spacing)
  cat("SIM decay: ", kern$decay, "\n", sep="")
  cat(spacing)
  cat("SIM Variance: ", kern$variance, "\n", sep="")
  cat(spacing)
  cat("RBF Variance: ", kern$rbf_variance, "\n", sep="")
  if(kern$gaussianInitial) {
    cat(spacing)
    cat("DISIM Initial Variance: ", kern$initialVariance, "\n", sep="")
  }
  ## cat(spacing)
  ## cat("SIM delay: %2.4f\n", kern$delay)
}


disimXdisimKernCompute <- function (disimKern1, disimKern2, t1, t2=t1) {
  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if ( disimKern1$inverseWidth != disimKern2$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern1$di_decay != disimKern2$di_decay)
    stop("Kernels cannot be cross combined if they have different driving input decays.")

  if ( disimKern1$di_variance != disimKern2$di_variance)
    stop("Kernels cannot be cross combined if they have different driving input variances.")

  if ( disimKern1$rbf_variance != disimKern2$rbf_variance)
    stop("Kernels cannot be cross combined if they have different RBF variances.")

  if ( disimKern1$isNormalised != disimKern2$isNormalised )
    stop("Both the DISIM kernels have to be either normalised or not.")

  l <- sqrt(2/disimKern1$inverseWidth)
  h1 <- .disimComputeH(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l)
  h2 <- .disimComputeH(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l)
  hp1 <- .disimComputeHPrime(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l)
  hp2 <- .disimComputeHPrime(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l)
  K <- h1 + t(h2) + hp1 + t(hp2)
  K <- 0.5*disimKern1$rbf_variance*disimKern1$di_variance*sqrt(disimKern1$variance)*sqrt(disimKern2$variance)*K

  if (!disimKern1$isNormalised)
    K <- K*l*sqrt(pi)

  if ("gaussianInitial" %in% names(disimKern1) && disimKern1$gaussianInitial && 
      "gaussianInitial" %in% names(disimKern2) && disimKern2$gaussianInitial) {
    if (disimKern1$initialVariance != disimKern2$initialVariance)
      stop("Kernels cannot be cross combined if they have different initial variances.");
    
    dim1 <- dim(as.matrix(t1))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t1, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    delta = disimKern1$di_decay;
    D1 = disimKern1$decay;
    D2 = disimKern2$decay;
  
    K = K + disimKern1$initialVariance * 
      sqrt(disimKern1$variance) * sqrt(disimKern2$variance) * 
        (exp(-delta * t1Mat) - exp(-D1 * t1Mat)) / (D1 - delta) *
          (exp(-delta * t2Mat) - exp(-D2 * t2Mat)) / (D2 - delta)
  }
  
  return (K)
}



disimXrbfKernCompute <- function (disimKern, rbfKern, t1, t2=t1) {

  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if ( disimKern$inverseWidth != rbfKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern$rbf_variance != rbfKern$variance )
    stop("Kernels cannot be cross combined if they have different RBF variances.")

  if ( disimKern$isNormalised != rbfKern$isNormalised )
    stop("Both kernels have to be either normalised or not.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]

  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  l <- sqrt(2/disimKern$inverseWidth)

  Di <- disimKern$decay
  delta <- disimKern$di_decay

  invLDiffT <- 1/l*diffT
  halfLD_i <- 0.5*l*Di
  halfLDelta <- 0.5*l*delta

  lnCommon <- - log(0i + delta - Di)
  lnPart1_f <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart2_f <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]
  lnPart2 <- lnPart2_f[[1]]
  signs2 <- lnPart2_f[[2]]

  K <- signs1 * exp(lnCommon + halfLDelta^2 - delta * diffT + lnPart1) + signs2 * exp(lnCommon + halfLD_i^2 - Di * diffT + lnPart2)

  K <- 0.5*sqrt(disimKern$variance)*sqrt(disimKern$di_variance)*rbfKern$variance*K
  if (!disimKern$isNormalised)
    K <- K*l*sqrt(pi)
  K <- Re(K)

  return (K)
}




disimXsimKernCompute <- function (disimKern, simKern, t1, t2=t1) {
  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if ( disimKern$inverseWidth != simKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern$di_decay != simKern$decay)
    stop("Kernels cannot be cross combined if they have different driving input decays.")

  if ( disimKern$di_variance != simKern$variance)
    stop("Kernels cannot be cross combined if they have different driving input variances.")

  if ( disimKern$isNormalised != simKern$isNormalised )
    stop("Both kernels have to be either normalised or not.")

#  if ( disimKern$rbf_variance != simKern$rbf_variance)
#    stop("Kernels cannot be cross combined if they have different RBF variances.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]

  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  l <- sqrt(2/disimKern$inverseWidth)

  Di <- disimKern$decay
  delta <- disimKern$di_decay

  invLDiffT <- 1/l*diffT
  halfLD_i <- 0.5*l*Di
  halfLDelta <- 0.5*l*delta

  lnCommon1 <- - log(2*delta) -delta * t2Mat - Di * t1Mat + halfLDelta^2

  lnFact1 <- log(2 * delta) - log(0i + delta^2 - Di^2)
  lnPart1_f <- lnDiffErfs(halfLDelta - t2Mat/l, halfLDelta)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]

  lnFact2 <- (Di - delta) * t1Mat - log(0i + delta - Di)
  lnPart2a_f <- lnDiffErfs(halfLDelta, halfLDelta - t1Mat/l)
  lnPart2b_f <- lnDiffErfs(halfLDelta, halfLDelta - t2Mat/l)

  lnPart2a <- lnPart2a_f[[1]]
  signs2a <- lnPart2a_f[[2]]
  lnPart2b <- lnPart2b_f[[1]]
  signs2b <- lnPart2b_f[[2]]

  lnFact3 <- (Di + delta) * t1Mat - log(delta + Di)
  lnPart3_f <- lnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT)
  lnPart3 <- lnPart3_f[[1]]
  signs3 <- lnPart3_f[[2]]

  lnFact4 <- 2*delta*t2Mat + (Di - delta) * t1Mat - log(0i + delta - Di)
  lnPart4_f <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart4 <- lnPart4_f[[1]]
  signs4 <- lnPart4_f[[2]]

  lnCommon2 <- - log(0i + delta^2 - Di^2) - delta * t2Mat - Di * t1Mat + halfLD_i^2
  lnPart5_f <- lnDiffErfs(halfLD_i - t1Mat/l, halfLD_i)
  lnPart5 <- lnPart5_f[[1]]
  signs5 <- lnPart5_f[[2]]

  lnFact6 <- (Di + delta) * t2Mat
  lnPart6_f <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)
  lnPart6 <- lnPart6_f[[1]]
  signs6 <- lnPart6_f[[2]]

  K <- signs1 * exp(lnCommon1 + lnFact1 + lnPart1) +signs2a*exp(lnCommon1 + lnFact2 + lnPart2a) +signs2b*exp(lnCommon1 + lnFact2 + lnPart2b) +signs3*exp(lnCommon1 + lnFact3 + lnPart3) +signs4*exp(lnCommon1 + lnFact4 + lnPart4) +signs5*exp(lnCommon2 + lnPart5) +signs6*exp(lnCommon2 + lnFact6 + lnPart6)
  K <- 0.5*disimKern$rbf_variance*disimKern$di_variance*sqrt(disimKern$variance)*K
  if (!disimKern$isNormalised)
    K <- K*l*sqrt(pi)
  K <- Re(K)

  if ("gaussianInitial" %in% names(disimKern) && disimKern$gaussianInitial && 
      "gaussianInitial" %in% names(simKern) && simKern$gaussianInitial) {
    if (disimKern$initialVariance != simKern$initialVariance)
      stop("Kernels cannot be cross combined if they have different initial variances.")
    
    dim1 <- dim(as.matrix(t1))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t1, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    delta = disimKern$di_decay
    D = disimKern$decay
  
    K = K + disimKern$initialVariance * 
      sqrt(disimKern$variance) *
        (exp(-delta * t1Mat) - exp(-D * t1Mat)) / (D - delta) *
          exp(-delta * t2Mat)
  }

  return (K)
}



disimKernDiagCompute <- function (kern, t1) {
  if ( dim(as.matrix(t1))[2]>1 )
    stop("Input can only have one column.")

  l <- sqrt(2/kern$inverseWidth)
  delta <- kern$di_decay
  D <- kern$decay
  halfLD <- 0.5*l*D
  halfLDelta <- 0.5*l*delta

  lnPart1_f <- lnDiffErfs(halfLDelta - t1/l, halfLDelta)
  lnPart2_f <- lnDiffErfs(halfLDelta + t1/l, halfLDelta)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]
  lnPart2 <- lnPart2_f[[1]]
  signs2 <- lnPart2_f[[2]]

  lnCommon <- halfLDelta ^ 2 -(D+delta)*t1 - log(2*delta) - log(0i + D-delta)
  lnFact2 <- (D+delta)*t1 - log(D + delta)

  if (abs(D - delta) < .1)
    h <- signs1 * exp(lnCommon + lnPart1) * ((exp((D-delta)*t1) - 1) / (D - delta) + 1/(D+delta)) + signs2 * exp(lnCommon + lnFact2 + lnPart2)
  else {
    lnFact1a <- (D - delta) * t1 + log(D + delta) - log(0i + D^2 - delta^2)
    lnFact1b <- log(2*delta) - log(0i + D^2 - delta^2)
    h <- signs1 * exp(lnCommon + lnFact1a + lnPart1) - signs1 * exp(lnCommon + lnFact1b + lnPart1) + signs2 * exp(lnCommon + lnFact2 + lnPart2)
  }

  lnPart1p_f <- lnDiffErfs(halfLD - t1/l, halfLD)
  lnPart2p_f <- lnDiffErfs(halfLD + t1/l, halfLD)
  lnPart1p <- lnPart1p_f[[1]]
  signs1p <- lnPart1p_f[[2]]
  lnPart2p <- lnPart2p_f[[1]]
  signs2p <- lnPart2p_f[[2]]

  lnCommonp <- halfLD^2 - 2*D*t1 - log(0i + delta^2 - D^2)
  lnFact2p <- 2*D*t1 - log(2*D)

  if (abs(D - delta) < .1) {
    hp <- signs1p * exp(lnCommonp + lnPart1p) * ((exp((D-delta)*t1) - 1) / (D - delta) + 1/(2*D)) + signs2p * exp(lnCommonp + lnFact2p + lnPart2p)
  }
  else {
    lnFact1ap <- log(D + delta) - log(0i + delta - D) - log(2*D)
    lnFact1bp <- (D-delta)*t1 - log(0i + delta - D)

    hp <- signs1p * exp(lnCommonp + lnFact1ap + lnPart1p) - signs1p * exp(lnCommonp + lnFact1bp + lnPart1p) + signs2p * exp(lnCommonp + lnFact2p + lnPart2p)
  }

  k <- kern$rbf_variance*kern$di_variance*kern$variance*Re(h+hp)
  if (!kern$isNormalised)
    k <- k*l*sqrt(pi)

  if ( "gaussianInitial" %in% names(kern) && kern$gaussianInitial )
    k = k + kern$initialVariance*kern$variance * 
      ((exp(-kern$di_decay*t1) - exp(-kern$decay*t1)) /
       (kern$decay-kern$di_decay))^2;

  return (k)
}



disimKernGradient <- function (kern, t1, t2, covGrad) {
  if ( nargs() == 3 ) {
    covGrad <- t2
    t2 <- t1
  }

  gFull <- disimXdisimKernGradient(kern, kern, t1, t2, covGrad)

  g <- gFull$g1 + gFull$g2

  return (g)
}



disimXdisimKernGradient <- function (disimKern1, disimKern2, t1, t2, covGrad) {
  if ( nargs() < 5 ) {
    covGrad <- t2
    t2 <- t1
  }

  if ( dim(as.matrix(t1))[2]>1 | dim(as.matrix(t2))[2]>1 )
    stop("Input can only have one column for SIM kernel.")

  if ( disimKern1$inverseWidth != disimKern2$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern1$di_decay != disimKern2$di_decay)
    stop("Kernels cannot be cross combined if they have different driving input decays.")

  if ( disimKern1$di_variance != disimKern2$di_variance )
    stop("Kernels cannot be cross combined if they have different driving input variances.")

  if ( disimKern1$rbf_variance != disimKern2$rbf_variance )
    stop("Kernels cannot be cross combined if they have different RBF variances.")

  if ( disimKern1$isNormalised != disimKern2$isNormalised )
    stop("Both the DISIM kernels have to be either normalised or not.")

  option <- 5
  l <- sqrt(2/disimKern1$inverseWidth)
  disimH1 <- .disimComputeH(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l, option)
  h1 <- disimH1$h
  dh1_ddelta <- disimH1$dh_ddelta
  dh1_dD1 <- disimH1$dh_dD_j
  dh1_dD2 <- disimH1$dh_dD_k
  dh1_dl <- disimH1$dh_dl

  disimH2 <- .disimComputeH(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l, option)
  h2 <- disimH2$h
  dh2_ddelta <- disimH2$dh_ddelta
  dh2_dD2 <- disimH2$dh_dD_j
  dh2_dD1 <- disimH2$dh_dD_k
  dh2_dl <- disimH2$dh_dl

  disimHp1 <- .disimComputeHPrime(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l, option)

  hp1 <- disimHp1$h
  dhp1_ddelta <- disimHp1$dh_ddelta
  dhp1_dD1 <- disimHp1$dh_dD_j
  dhp1_dD2 <- disimHp1$dh_dD_k
  dhp1_dl <- disimHp1$dh_dl

  disimHp2 <- .disimComputeHPrime(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l, option)
  hp2 <- disimHp2$h
  dhp2_ddelta <- disimHp2$dh_ddelta
  dhp2_dD2 <- disimHp2$dh_dD_j
  dhp2_dD1 <- disimHp2$dh_dD_k
  dhp2_dl <- disimHp2$dh_dl

  dK_ddelta <- dh1_ddelta + t(dh2_ddelta) + dhp1_ddelta + t(dhp2_ddelta)
  dK_dD1 <- dh1_dD1 + t(dh2_dD1) + dhp1_dD1 + t(dhp2_dD1)
  dK_dD2 <- dh1_dD2 + t(dh2_dD2) + dhp1_dD2 + t(dhp2_dD2)
  dK_dl <- dh1_dl + t(dh2_dl) + dhp1_dl + t(dhp2_dl)

  C0 <- disimKern1$di_variance
  C1 <- sqrt(disimKern1$variance)
  C2 <- sqrt(disimKern2$variance)
  C3 <- disimKern1$rbf_variance
  K <- 0.5 * (h1 + t(h2) + hp1 + t(hp2))
  var2 <- C0*C1*C2*C3

  if (disimKern1$isNormalised) {
    dk_ddelta <- sum(covGrad*dK_ddelta)*0.5*var2
    dk_dD1 <- sum(covGrad*dK_dD1)*0.5*var2
    dk_dD2 <- sum(covGrad*dK_dD2)*0.5*var2
    dk_dl <- sum(covGrad*dK_dl)*0.5*var2
  }
  else {
    K <- K*sqrt(pi)
    dk_ddelta <- sum(covGrad*dK_ddelta)*0.5*sqrt(pi)*l*var2
    dk_dD1 <- sum(covGrad*dK_dD1)*0.5*sqrt(pi)*l*var2
    dk_dD2 <- sum(covGrad*dK_dD2)*0.5*sqrt(pi)*l*var2
    dk_dl <- sum(covGrad*(dK_dl*0.5*sqrt(pi)*l + K))*var2
    K <- l*K
  }
  dk_dC0 <- C1*C2*C3*sum(covGrad*K)
  dk_dC1 <- C0*C2*C3*sum(covGrad*K)
  dk_dC2 <- C0*C1*C3*sum(covGrad*K)
  dk_dC3 <- C0*C1*C2*sum(covGrad*K)

  dk_dDIVariance <- dk_dC0
  dk_dDisim1Variance <- dk_dC1*0.5/C1
  dk_dDisim2Variance <- dk_dC2*0.5/C2
  dk_dRBFVariance <- dk_dC3

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern1$inverseWidth*sqrt(disimKern1$inverseWidth))*dk_dl

  K <- var2*K

  if ("gaussianInitial" %in% names(disimKern1) && disimKern1$gaussianInitial && 
      "gaussianInitial" %in% names(disimKern2) && disimKern2$gaussianInitial) {
    if (disimKern1$initialVariance != disimKern2$initialVariance)
      stop("Kernels cannot be cross combined if they have different initial variances.");
    
    dim1 <- dim(as.matrix(t1))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t1, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    delta <- disimKern1$di_decay
    D1 <- disimKern1$decay
    D2 <- disimKern2$decay

    the_rest <- (exp(- delta * t1Mat) - exp(- D1 * t1Mat)) / (D1 - delta) *
      (exp(- delta * t2Mat) - exp(- D2 * t2Mat)) / (D2 - delta)
  
    dk_dinitVariance <- sum((sqrt(disimKern1$variance) *
      sqrt(disimKern2$variance) * the_rest) * covGrad)

    dk_dDisim1Variance <- dk_dDisim1Variance + 
      sum((.5 / sqrt(disimKern1$variance) *
           disimKern1$initialVariance * sqrt(disimKern2$variance) * the_rest)
          * covGrad)

    dk_dDisim2Variance <- dk_dDisim2Variance +
      sum((.5 / sqrt(disimKern2$variance) *
           disimKern1$initialVariance * sqrt(disimKern1$variance) * the_rest)
          * covGrad)

    dk_dD1 <- dk_dD1 +
      sum((disimKern1$initialVariance *
           sqrt(disimKern1$variance) * sqrt(disimKern2$variance) *
           (t1Mat * (D1 - delta)*exp(-D1*t1Mat) - exp(-delta*t1Mat) + 
            exp(-D1*t1Mat)) / (D1-delta)^2 *
           (exp(- delta * t2Mat) - exp(- D2 * t2Mat)) / (D2 - delta))
          *covGrad)
  
    dk_dD2 <- dk_dD2 +
      sum((disimKern1$initialVariance *
           sqrt(disimKern1$variance) * sqrt(disimKern2$variance) *
           (t2Mat * (D2 - delta)*exp(-D2*t2Mat) - exp(-delta*t2Mat) + 
            exp(-D2*t2Mat)) / (D2-delta)^2 *
           (exp(- delta * t1Mat) - exp(-D1 * t1Mat)) / (D1 - delta))
          *covGrad)

    dk_ddelta <- dk_ddelta +
      sum((disimKern1$initialVariance *
           sqrt(disimKern1$variance) * sqrt(disimKern2$variance) *
           ((-t2Mat * (D2 - delta)*exp(-delta*t2Mat) + exp(-delta*t2Mat) - 
             exp(-D2*t2Mat)) / (D2-delta)^2 *
            (exp(- delta * t1Mat) - exp(-D1 * t1Mat)) / (D1 - delta) +
            (-t1Mat * (D1 - delta)*exp(-delta*t1Mat) + exp(-delta*t1Mat) - 
             exp(-D1*t1Mat)) / (D1-delta)^2 *
            (exp(- delta * t2Mat) - exp(-D2 * t2Mat)) / (D2 - delta)))
          *covGrad)
    
    g1 <- c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD1, dk_dDisim1Variance, dk_dRBFVariance, dk <- dk_dinitVariance)
    g2 <- c(0, 0, 0, dk_dD2, dk_dDisim2Variance, 0, 0)

  }
  else {
    g1 <- c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD1, dk_dDisim1Variance, dk_dRBFVariance)
    g2 <- c(0, 0, 0, dk_dD2, dk_dDisim2Variance, 0)
  }

  g <- list(g1=g1, g2=g2)
  return (g)
}



disimXrbfKernGradient <- function (disimKern, rbfKern, t1, t2, covGrad) {
  if ( nargs() < 5 ) {
    covGrad <- t2
    t2 <- t1
  }

  if ( dim(as.matrix(t1))[2]>1 | dim(as.matrix(t2))[2]>1 )
    stop("Input can only have one column for SIM kernel.")

  if ( disimKern$inverseWidth != rbfKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern$rbf_variance != rbfKern$variance )
    stop("Kernels cannot be cross combined if they have different RBF variances.")

  if ( disimKern$isNormalised != rbfKern$isNormalised )
    stop("Both kernels have to be either normalised or not.")

  if ( nargs()<5 ) {
    k <- disimXrbfKernCompute(disimKern, rbfKern, t1)
  } else {
    k <- disimXrbfKernCompute(disimKern, rbfKern, t1, t2)
  }

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  l <- sqrt(2/disimKern$inverseWidth)
  l2 <- l*l
  C_0 <- sqrt(disimKern$di_variance)
  C_i <- sqrt(disimKern$variance)
  C_j <- rbfKern$variance
  D_i <- disimKern$decay
  delta <-  disimKern$di_decay

  invLDiffT <- 1/l*diffT
  halfLD_i <- 0.5*l*D_i
  halfLDelta <- 0.5*l*delta

  if (disimKern$isNormalised)
    prefact <- 0.5 * C_0 * C_i * C_j
  else
    prefact <- C_0 * C_i * C_j * sqrt(pi)/2 * l

  lnCommon <- - log(0i + delta - D_i)
  lnPart1_f <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart2_f <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)
  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]
  lnPart2 <- lnPart2_f[[1]]
  signs2 <- lnPart2_f[[2]]

  lnFact1 <- halfLDelta^2 - delta * diffT
  lnFact2 <- halfLD_i^2 - D_i * diffT

  gradln1 <- .gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, l/2, l/2)
  dlnPart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2 <-.gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, l/2, l/2)
  dlnPart2 <- gradln2$dlnPart
  m2 <- gradln2$m

  dK_dD <- k * (1/(delta-D_i)) + prefact * ((l*halfLD_i - diffT) * signs2 * exp(lnCommon + lnFact2 + lnPart2) + dlnPart2 * exp(lnCommon + lnFact2 - m2))
  dk_dD <- sum(dK_dD*covGrad)

  dK_ddelta <- k * (-1/(delta-D_i)) + prefact * ((l*halfLDelta - diffT) * signs1 * exp(lnCommon + lnFact1 + lnPart1) + dlnPart1 * exp(lnCommon + lnFact1 - m1))
  dk_ddelta <- sum(dK_ddelta*covGrad)

  gradln1 <- .gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, delta/2 + invLDiffT/l, delta/2 - t2Mat/l2)
  dlnPart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2 <-.gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, D_i/2 - t2Mat/l2, D_i/2 + invLDiffT/l)
  dlnPart2 <- gradln2$dlnPart
  m2 <- gradln2$m

  dK_dl <- prefact * (delta*halfLDelta * signs1 * exp(lnCommon + lnFact1 + lnPart1) + D_i*halfLD_i * signs2 * exp(lnCommon + lnFact2 + lnPart2) + dlnPart1 * exp(lnCommon + lnFact1 - m1) + dlnPart2 * exp(lnCommon + lnFact2 - m2))
  if (!disimKern$isNormalised)
    dK_dl <- dK_dl + k/l
  
  dk_dl <- sum(dK_dl*covGrad)

  dk_dC_i <- sum(k*covGrad)/C_i
  dk_dC_0 <- sum(k*covGrad)/C_0
  dk_dRbfVariance <- sum(k*covGrad)/rbfKern$variance

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern$inverseWidth*sqrt(disimKern$inverseWidth))*dk_dl
  dk_dDisimVariance <- dk_dC_i*0.5/C_i
  dk_dDIVariance <- dk_dC_0*0.5/C_0

  if ("gaussianInitial" %in% names(disimKern) && disimKern$gaussianInitial)
    g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dDisimVariance, 0, 0))
  else
    g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dDisimVariance, 0))
  g2 <- Re(c(0, dk_dRbfVariance))

  g <- list(g1=g1, g2=g2)
  return (g)
}




disimXsimKernGradient <- function (disimKern, simKern, t1, t2, covGrad) {
  if ( nargs() < 5 ) {
    covGrad <- t2
    t2 <- t1
  }

  if ( dim(as.matrix(t1))[2]>1 | dim(as.matrix(t2))[2]>1 )
    stop("Input can only have one column for SIM kernel.")

  if ( disimKern$inverseWidth != simKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

  if ( disimKern$di_decay != simKern$decay)
    stop("Kernels cannot be cross combined if they have different driving input decays.")

  if ( disimKern$di_variance != simKern$variance)
    stop("Kernels cannot be cross combined if they have different driving input variances.")

  if ( disimKern$isNormalised != simKern$isNormalised )
    stop("Both kernels have to be either normalised or not.")

#  if ( disimKern$rbf_variance != simKern$rbf_variance)
#    stop("Kernels cannot be cross combined if they have different RBF variances.")

  dim1 <- dim(as.matrix(t1))[1]
  dim2 <- dim(as.matrix(t2))[1]
  t1Mat <- matrix(t1, dim1, dim2)
  t2Mat <- t(matrix(t2, dim2, dim1))
  diffT <- t1Mat-t2Mat

  l <- sqrt(2/disimKern$inverseWidth)
  l2 <- l*l
  C_0 <- sqrt(disimKern$di_variance)
  C_i <- sqrt(disimKern$variance)

  D_i <- disimKern$decay
  delta <-  disimKern$di_decay
  halfLD_i <- 0.5*l*D_i
  halfLDelta <- 0.5*l*delta
  invLDiffT <- 1/l*diffT

  lnCommon1 <- - log(2*delta) -delta * t2Mat - D_i * t1Mat + halfLDelta^2

  lnFact1 <- log(2 * delta) - log(0i + delta^2 - D_i^2)
  lnPart1_f <- lnDiffErfs(halfLDelta - t2Mat/l, halfLDelta)

  lnPart1 <- lnPart1_f[[1]]
  signs1 <- lnPart1_f[[2]]

  lnFact2 <- (D_i - delta) * t1Mat - log(0i + delta - D_i)
  lnPart2a_f <- lnDiffErfs(halfLDelta, halfLDelta - t1Mat/l)
  lnPart2b_f <- lnDiffErfs(halfLDelta, halfLDelta - t2Mat/l)

  lnPart2a <- lnPart2a_f[[1]]
  signs2a <- lnPart2a_f[[2]]
  lnPart2b <- lnPart2b_f[[1]]
  signs2b <- lnPart2b_f[[2]]

  lnFact3 <- (D_i + delta) * t1Mat - log(delta + D_i)
  lnPart3_f <- lnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT)
  lnPart3 <- lnPart3_f[[1]]
  signs3 <- lnPart3_f[[2]]

  lnFact4 <- 2*delta*t2Mat + (D_i - delta) * t1Mat - log(0i + delta - D_i)
  lnPart4_f <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart4 <- lnPart4_f[[1]]
  signs4 <- lnPart4_f[[2]]

  lnCommon2 <- - log(0i + delta^2 - D_i^2) - delta * t2Mat - D_i * t1Mat + halfLD_i^2
  lnPart5_f <- lnDiffErfs(halfLD_i - t1Mat/l, halfLD_i)
  lnPart5 <- lnPart5_f[[1]]
  signs5 <- lnPart5_f[[2]]

  lnFact6 <- (D_i + delta) * t2Mat
  lnPart6_f <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)
  lnPart6 <- lnPart6_f[[1]]
  signs6 <- lnPart6_f[[2]]

  K1 <- signs1 * exp( lnCommon1 + lnFact1 + lnPart1) +signs2a*exp(lnCommon1 + lnFact2 + lnPart2a) +signs2b*exp(lnCommon1 + lnFact2 + lnPart2b) +signs3*exp(lnCommon1 + lnFact3 + lnPart3) +signs4*exp(lnCommon1 + lnFact4 + lnPart4)

  K2 <- signs5*exp( lnCommon2 + lnPart5) +signs6*exp(lnCommon2 + lnFact6 + lnPart6)

  if (disimKern$isNormalised)
    prefact <- 0.5*disimKern$rbf_variance*disimKern$di_variance*sqrt(disimKern$variance)
  else
    prefact <- 0.5*sqrt(pi)*l*disimKern$rbf_variance*disimKern$di_variance*sqrt(disimKern$variance)
  K <- prefact*Re(K1+K2)

  dcommon1 <- - 1/delta - t2Mat + l*halfLDelta
  dfact1 <- 1/delta - 2*delta/(delta^2 - D_i^2)
  dfact2 <- -t1Mat - 1/(delta - D_i)
  dfact3 <- t1Mat - 1/(delta + D_i)
  dfact4 <- 2*t2Mat - t1Mat - 1/(delta - D_i)
  dcommon2 <- - 2*delta/(delta^2 - D_i^2) - t2Mat
  dfact6 <- t2Mat

  gradln1 <- .gradLnDiffErfs(halfLDelta - t2Mat/l, halfLDelta, l/2, l/2)
  dpart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2a <-.gradLnDiffErfs(halfLDelta, halfLDelta - t1Mat/l, l/2, l/2)
  dpart2a <- gradln2a$dlnPart
  m2a <- gradln2a$m

  gradln2b <-.gradLnDiffErfs(halfLDelta, halfLDelta - t2Mat/l, l/2, l/2)
  dpart2b <- gradln2b$dlnPart
  m2b <- gradln2b$m

  gradln3 <- .gradLnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT, l/2, l/2)
  dpart3 <- gradln3$dlnPart
  m3 <- gradln3$m

  gradln4 <- .gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, l/2, l/2)
  dpart4 <- gradln4$dlnPart
  m4 <- gradln4$m

  dK_ddelta <- prefact * (dcommon1 * K1 + dcommon2 * K2 + dfact1 * signs1*exp(lnCommon1 + lnFact1 + lnPart1) + dfact2 * signs2a*exp(lnCommon1 + lnFact2 + lnPart2a) + dfact2 * signs2b*exp(lnCommon1 + lnFact2 + lnPart2b) + dfact3 * signs3*exp(lnCommon1 + lnFact3 + lnPart3) + dfact4 * signs4*exp(lnCommon1 + lnFact4 + lnPart4) + dfact6 * signs6*exp(lnCommon2 + lnFact6 + lnPart6) + dpart1 * exp(lnCommon1 + lnFact1 - m1) + dpart2a * exp(lnCommon1 + lnFact2 - m2a) + dpart2b * exp(lnCommon1 + lnFact2 - m2b) + dpart3 * exp(lnCommon1 + lnFact3 - m3) + dpart4 * exp(lnCommon1 + lnFact4 - m4))

  dcommon1 <- delta * halfLDelta
  dcommon2 <- D_i * halfLD_i

  gradln1 <- .gradLnDiffErfs(halfLDelta - t2Mat/l, halfLDelta, delta/2 + t2Mat/l2, delta/2)
  dpart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2a <-.gradLnDiffErfs(halfLDelta, halfLDelta - t1Mat/l, delta/2, delta/2 + t1Mat/l2)
  dpart2a <- gradln2a$dlnPart
  m2a <- gradln2a$m

  gradln2b <-.gradLnDiffErfs(halfLDelta, halfLDelta - t2Mat/l, delta/2, delta/2 + t2Mat/l2)
  dpart2b <- gradln2b$dlnPart
  m2b <- gradln2b$m

  gradln3 <- .gradLnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT, delta/2 - t1Mat/l2, delta/2 - invLDiffT/l)
  dpart3 <- gradln3$dlnPart
  m3 <- gradln3$m

  gradln4 <- .gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, delta/2 + invLDiffT/l, delta/2 - t2Mat/l2)
  dpart4 <- gradln4$dlnPart
  m4 <- gradln4$m

  gradln5 <- .gradLnDiffErfs(halfLD_i - t1Mat/l, halfLD_i, D_i/2 + t1Mat/l2, D_i/2)
  dpart5 <- gradln5$dlnPart
  m5 <- gradln5$m

  gradln6 <- .gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, D_i/2 - t2Mat/l2, D_i/2 + invLDiffT/l)
  dpart6 <- gradln6$dlnPart
  m6 <- gradln6$m

  dK_dl <- Re(prefact * (dcommon1 * K1 + dcommon2 * K2 +dpart1 * exp(lnCommon1 + lnFact1 - m1) +dpart2a * exp(lnCommon1 + lnFact2 - m2a) +dpart2b * exp(lnCommon1 + lnFact2 - m2b) +dpart3 * exp(lnCommon1 + lnFact3 - m3) +dpart4 * exp(lnCommon1 + lnFact4 - m4) +dpart5 * exp(lnCommon2 - m5) +dpart6 * exp(lnCommon2 + lnFact6 - m6)))
  if (!disimKern$isNormalised)
    dK_dl <- dK_dl + K / l

  dcommon1 <- - t1Mat
  dfact1 <- 2 * D_i / (delta^2 - D_i^2)
  dfact2 <- t1Mat + 1/(delta - D_i)
  dfact3 <- t1Mat - 1/(delta + D_i)
  dfact4 <- t1Mat + 1/(delta - D_i)
  dcommon2 <- 2 * D_i / (delta^2 - D_i^2) - t1Mat + l*halfLD_i
  dfact6 <- t2Mat

  gradln5 <- .gradLnDiffErfs(halfLD_i - t1Mat/l, halfLD_i, l/2, l/2)
  dpart5 <- gradln5$dlnPart
  m5 <- gradln5$m

  gradln6 <- .gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, l/2, l/2)
  dpart6 <- gradln6$dlnPart
  m6 <- gradln6$m

  dK_dD = prefact * (dcommon1 * K1 + dcommon2 * K2 + dfact1 * signs1*exp(lnCommon1 + lnFact1 + lnPart1) + dfact2 * signs2a*exp(lnCommon1 + lnFact2 + lnPart2a) + dfact2 * signs2b*exp(lnCommon1 + lnFact2 + lnPart2b) + dfact3 * signs3*exp(lnCommon1 + lnFact3 + lnPart3) + dfact4 * signs4*exp(lnCommon1 + lnFact4 + lnPart4) + dfact6 * signs6*exp(lnCommon2 + lnFact6 + lnPart6) + dpart5 * exp(lnCommon2 - m5) + dpart6 * exp(lnCommon2 + lnFact6 - m6))

  dk_ddelta <- sum(dK_ddelta*covGrad)
  dk_dl <- sum(dK_dl*covGrad)
  dk_dD <- sum(dK_dD*covGrad)

  dk_dRBFVariance <- sum(K*covGrad)/disimKern$rbf_variance
  dk_dDIVariance <- sum(K*covGrad)/disimKern$di_variance
  dk_dSimVariance <- .5 * sum(K*covGrad)/disimKern$variance

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern$inverseWidth* sqrt(disimKern$inverseWidth))*dk_dl

  if ("gaussianInitial" %in% names(disimKern) && disimKern$gaussianInitial && 
      "gaussianInitial" %in% names(simKern) && simKern$gaussianInitial) {
    if (disimKern$initialVariance != simKern$initialVariance)
      stop("Kernels cannot be cross combined if they have different initial variances.")
    
    dim1 <- dim(as.matrix(t1))[1]
    dim2 <- dim(as.matrix(t2))[1]
    t1Mat <- matrix(t1, dim1, dim2)
    t2Mat <- t(matrix(t2, dim2, dim1))

    delta = disimKern$di_decay
    D = disimKern$decay

    the_rest = (exp(-delta*t1Mat) - exp(-D*t1Mat)) / (D-delta) *
      exp(-delta*t2Mat)
  
    dk_dinitVariance = sum((sqrt(disimKern$variance) * the_rest) * covGrad)

    dk_dSimVariance = dk_dSimVariance +
      sum((.5 / sqrt(disimKern$variance) *
           disimKern$initialVariance * the_rest) * covGrad)

    dk_dD = dk_dD +
      sum((disimKern$initialVariance * sqrt(disimKern$variance) *
           (t1Mat*(D-delta)*exp(-D*t1Mat) - exp(-delta*t1Mat)
            + exp(-D*t1Mat)) / (D-delta)^2 *
           exp(-delta*t2Mat))*covGrad)
  
    dk_ddelta = dk_ddelta +
      sum((disimKern$initialVariance * sqrt(disimKern$variance) *
           (-t2Mat*exp(-delta*t2Mat) *
            (exp(-delta*t1Mat) - exp(-D*t1Mat)) / (D-delta) +
            (-t1Mat*(D-delta)*exp(-delta*t1Mat) + exp(-delta*t1Mat)
             - exp(-D*t1Mat)) / (D-delta)^2 *
            exp(-delta*t2Mat)))*covGrad)
    
    g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dSimVariance, dk_dRBFVariance, dk_dinitVariance))
    g2 <- c(0, 0, 0, 0)
  }
  else {
    g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dSimVariance, dk_dRBFVariance))
    g2 <- c(0, 0, 0)
  }

  g <- list(g1=g1, g2=g2)
  return (g)
}
