disimKernParamInit <- function (kern) {

  if ( kern$inputDimension > 1 )
    stop("SIM kernel is only valid for one-D input.")

  kern$di_decay <- 0.1
  kern$inverseWidth <- 1
  kern$di_variance <- 1
  kern$decay <- 1
  kern$variance <- 1
  kern$nParams <- 5

  kern$transforms <- list(index=c(1,2,3,4,5), type="positive")

  kern$isStationary=FALSE

  return (kern)

}



disimKernCompute <- function (kern, t, t2=t) {

  if ( ( dim(as.matrix(t))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  l <- sqrt(2/kern$inverseWidth)

  h <- disimComputeH(t, t2, kern$di_decay, kern$decay, kern$decay, l)
  hp <- disimComputeHPrime(t, t2, kern$di_decay, kern$decay, kern$decay, l)
  if ( nargs()<3 ) {
    k <- h + t(h) + hp + t(hp)
  } else {
    h2 <- disimComputeH(t2, t, kern$di_decay, kern$decay, kern$decay, l)
    hp2 <- disimComputeHPrime(t2, t, kern$di_decay, kern$decay, kern$decay, l)
    k <- h + t(h2) + hp + t(hp2)
  }
  k <- 0.5*k*sqrt(pi)*l
  k <- kern$di_variance*kern$variance*k
  k <- Re(k)

  return (k)
}



disimComputeH <-  function (t1, t2, delta, Dj, Dk, l, option=1) {

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

  lnPart1 <- lnDiffErfs(halfLDelta - t1Mat/l, halfLDelta)
  lnPart2 <- lnDiffErfs(halfLDelta + t2Mat/l, halfLDelta-invLDiffT)

  lnCommon <- halfLDelta^2-complexLog(2*delta)-Dk*t2Mat-delta*t1Mat-complexLog(Dj-delta)

  lnFact1a1 <- complexLog(Dk + delta) + (Dk-delta)*t2Mat
  lnFact1a2 <- complexLog(2*delta)
  lnFact1b <- complexLog(Dk^2 - delta^2)

  lnFact2a <- (Dk + delta)*t2Mat
  lnFact2b <- complexLog(Dk + delta)

  h <- exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) - exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) + exp(lnCommon + lnFact2a - lnFact2b + lnPart2)

  h <- Re(h)

  l2 <- l*l

  if ( option == 1 ) {
    return (h)
  } else {
    m1 <- apply((halfLDelta - t1Mat/l)^2, c(1,2), min, halfLDelta^2)
    m2 <- apply((halfLDelta + t2Mat/l)^2, c(1,2), min, (halfLDelta - invLDiffT)^2)

    dlnPart1 <- l/sqrt(pi)*(exp(-(halfLDelta - t1Mat/l)^2 + m1) - exp(-halfLDelta^2 + m1))
    dlnPart2 <- l/sqrt(pi)*(exp(-(halfLDelta + t2Mat/l)^2 + m2) - exp(-(halfLDelta - invLDiffT)^2 + m2))
    dh_ddelta <- (l*halfLDelta - t1Mat + 1/(Dj-delta)-1/delta)*h + exp(lnCommon + lnFact1a1 - lnFact1b - m1)*dlnPart1 - exp(lnCommon + lnFact1a2 - lnFact1b - m1)*dlnPart1 + exp(lnCommon + lnFact2a - lnFact2b - m2)*dlnPart2 + exp(lnCommon + lnPart1 + lnFact1a1 - lnFact1b)*(1/(Dk-delta) - t2Mat) - exp(lnCommon + lnPart1 + complexLog(2*(Dk^2 + delta^2)) -2*lnFact1b) + (t2Mat - 1/(Dk + delta))*exp((Dk + delta)*t2Mat + lnCommon + lnPart2 - complexLog(Dk + delta))
    dh_ddelta <- Re(dh_ddelta)

    disimH <- list(h=h, dh_ddelta=dh_ddelta)

    if ( option > 2 ) {
      dh_dD_j = Re(- 1/(Dj - delta)*h)

      disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j)

      if ( option > 3 ) {
        dh_dD_k <- -t2Mat*h + exp(lnCommon + lnPart1 + lnFact1a1 - lnFact1b) * (t2Mat - 1/(Dk - delta)) - exp(lnCommon + lnPart1) * (-4*delta*Dk / (Dk^2 - delta^2)^2) + (t2Mat - 1/(Dk + delta)) * exp(lnFact2a - lnFact2b + lnCommon + lnPart2)
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



disimComputeHPrime <-  function (t1, t2, delta, Dj, Dk, l, option=1) {

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

  lnPart1 <- lnDiffErfs(halfLD_k - t2Mat/l, halfLD_k)
  lnPart2 <- lnDiffErfs(halfLD_k + t1Mat/l, halfLD_k+invLDiffT)

  lnCommon <- halfLD_k^2 - Dj*t1Mat - Dk*t2Mat - complexLog(delta^2-Dk^2) - complexLog(Dj+Dk)

  lnFact1a1 <- complexLog(Dk + delta)
  lnFact1a2 <- complexLog(Dj + Dk) + (Dj-delta)*t1Mat
  lnFact1b <- complexLog(delta-Dj)

  lnFact2a <- (Dj + Dk) * t1Mat

  h <- exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) - exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) + exp(lnCommon + lnFact2a + lnPart2)

  h <- Re(h)

  l2 <- l*l

  if ( option == 1 ) {
    return (h)
  } else {
    dh_ddelta <- -2*delta/(delta^2-Dk^2)*h - exp(lnCommon + lnPart1 + complexLog(Dk + Dj) - 2*lnFact1b) - exp(lnCommon + lnPart1 + lnFact1a2 - lnFact1b)*(-t1Mat - 1/(delta - Dj))
    dh_ddelta <- Re(dh_ddelta)
    disimH <- list(h=h, dh_ddelta=dh_ddelta)

    if ( option > 2 ) {
      dh_dD_j <- (-t1Mat - 1 /(Dj + Dk))*h + exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1)*(1/(delta - Dj)) - exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1)*(1/(Dj + Dk) + t1Mat + 1/(delta - Dj)) + t1Mat*exp(lnCommon + lnPart2 + lnFact2a)
      dh_dD_j <- Re(dh_dD_j)

      disimH <- list(h=h, dh_ddelta=dh_ddelta, dh_dD_j=dh_dD_j)
      if ( option > 3 ) {
        m1 <- apply((halfLD_k - t2Mat/l)^2, c(1,2), min, halfLD_k^2)
        m2 <- apply((halfLD_k + t1Mat/l)^2, c(1,2), min, (halfLD_k - invLDiffT)^2)

        dlnPart1 <- l/sqrt(pi) * (exp(-(halfLD_k - t2Mat/l)^2 + m1) - exp(-halfLD_k^2 + m1))
        dlnPart2 <- l/sqrt(pi) * (exp(-(halfLD_k + t1Mat/l)^2 + m2) - exp(-(halfLD_k + invLDiffT)^2 + m2))
        dh_dD_k <- (l*halfLD_k + 2*Dk / (delta^2 - Dk^2) + (-t2Mat - 1/(Dj + Dk)))*h + exp(lnCommon + lnFact1a1 - lnFact1b - m1) * dlnPart1 - exp(lnCommon + lnFact1a2 - lnFact1b - m1) * dlnPart1 + exp(lnCommon + lnFact2a - m2) * dlnPart2 + exp(lnCommon + lnFact1a1 - lnFact1b + lnPart1) * (1/(Dk + delta)) - exp(lnCommon + lnFact1a2 - lnFact1b + lnPart1) *(1/(Dj + Dk)) + exp(lnCommon + lnPart2 + complexLog(t1Mat) + lnFact2a)
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



disimKernExtractParam <- function (kern, option=1) {

  if ( option == 1 ) {
    params <- c(kern$di_decay, kern$inverseWidth, kern$di_variance, kern$decay, kern$variance)

  } else {
    params <- list(values=c(kern$di_decay, kern$inverseWidth, kern$di_variance, kern$decay, kern$variance), names=c("di_decay", "inverseWidth", "di_variance", "decay", "variance"))
  }

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

  return (kern)
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

  l <- sqrt(2/disimKern1$inverseWidth)
  h1 <- disimComputeH(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l)
  h2 <- disimComputeH(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l)
  hp1 <- disimComputeHPrime(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l)
  hp2 <- disimComputeHPrime(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l)
  K <- h1 + t(h2) + hp1 + t(hp2)
  K <- 0.5*K*sqrt(pi)*l
  K <- disimKern1$di_variance*sqrt(disimKern1$variance)*sqrt(disimKern2$variance)*K

  return (K)
}



disimXrbfKernCompute <- function (disimKern, rbfKern, t1, t2=t1) {

  if ( ( dim(as.matrix(t1))[2] > 1 ) | ( dim(as.matrix(t2))[2] > 1 ) )
    stop("Input can only have one column.")

  if ( disimKern$inverseWidth != rbfKern$inverseWidth )
    stop("Kernels cannot be cross combined if they have different inverse widths.")

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

  lnCommon <- - complexLog(delta - Di)
  lnPart1 <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart2 <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)

  K <- exp(lnCommon + halfLDelta^2 - delta * diffT + lnPart1) + exp(lnCommon + halfLD_i^2 - Di * diffT + lnPart2)

  K <- 0.5*sqrt(disimKern$variance)*sqrt(disimKern$di_variance)*sqrt(rbfKern$variance)*K*sqrt(pi)*l
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

  lnCommon1 <- - complexLog(2*delta) -delta * t2Mat - Di * t1Mat + halfLDelta^2

  lnFact1 <- complexLog(2 * delta) - complexLog(delta^2 - Di^2)
  lnPart1 <- lnDiffErfs(halfLDelta - t2Mat/l, halfLDelta)

  lnFact2 <- (Di - delta) * t1Mat - complexLog(delta - Di)
  lnPart2a <- lnDiffErfs(halfLDelta, halfLDelta - t1Mat/l)
  lnPart2b <- lnDiffErfs(halfLDelta, halfLDelta - t2Mat/l)

  lnFact3 <- (Di + delta) * t1Mat - complexLog(delta + Di)
  lnPart3 <- lnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT)

  lnFact4 <- 2*delta*t2Mat + (Di - delta) * t1Mat - complexLog(delta - Di)
  lnPart4 <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)

  lnCommon2 <- - complexLog(delta^2 - Di^2) - delta * t2Mat - Di * t1Mat + halfLD_i^2
  lnPart5 <- lnDiffErfs(halfLD_i - t1Mat/l, halfLD_i)

  lnFact6 <- (Di + delta) * t2Mat
  lnPart6 <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)

  K <- exp(lnCommon1 + lnFact1 + lnPart1) +exp(lnCommon1 + lnFact2 + lnPart2a) +exp(lnCommon1 + lnFact2 + lnPart2b) +exp(lnCommon1 + lnFact3 + lnPart3) +exp(lnCommon1 + lnFact4 + lnPart4) +exp(lnCommon2 + lnPart5) +exp(lnCommon2 + lnFact6 + lnPart6)
  K <- 0.5*K*sqrt(pi)*l
  K <- disimKern$di_variance*sqrt(disimKern$variance)*K
  K <- Re(K)

  return (K)
}



disimKernDiagCompute <- function (kern, t) {
  if ( dim(as.matrix(t))[2]>1 )
    stop("Input can only have one column.")

  l <- sqrt(2/kern$inverseWidth)
  delta <- kern$di_decay
  D <- kern$decay
  halfLD <- 0.5*l*D
  halfLDelta <- 0.5*l*delta

  lnPart1 <- lnDiffErfs(halfLDelta - t/l, halfLDelta)
  lnPart2 <- lnDiffErfs(halfLDelta + t/l, halfLDelta)

  lnCommon <- halfLDelta ^ 2 -(D+delta)*t - complexLog(2*delta) - complexLog(D-delta)
  lnFact1a <- (D - delta) * t + complexLog(D + delta) - complexLog(D^2 - delta^2)
  lnFact1b <- complexLog(2*delta) - complexLog(D^2 - delta^2)
  lnFact2 <- (D+delta)*t - complexLog(D + delta)

  h <- exp(lnCommon + lnFact1a + lnPart1) - exp(lnCommon + lnFact1b + lnPart1) + exp(lnCommon + lnFact2 + lnPart2)

  lnPart1p <- lnDiffErfs(halfLD - t/l, halfLD)
  lnPart2p <- lnDiffErfs(halfLD + t/l, halfLD)

  lnCommonp <- halfLD^2 - 2*D*t - complexLog(2*D) - complexLog(delta^2 - D^2)
  lnFact1ap <- complexLog(D + delta) - complexLog(delta - D)
  lnFact1bp <- complexLog(2*D) + (D-delta)*t - complexLog(delta - D)
  lnFact2p <- 2*D*t

  hp <- exp(lnCommonp + lnFact1ap + lnPart1p) - exp(lnCommonp + lnFact1bp + lnPart1p) + exp(lnCommonp + lnFact2p + lnPart2p)

  k <- 2*Re(h+hp)
  k <- 0.5*k*sqrt(pi)*l
  k <- kern$di_variance*kern$variance*k

  return (k)
}



disimKernGradient <- function (kern, t, t2, covGrad) {
  if ( nargs() == 3 ) {
    covGrad <- t2
    t2 <- t
  }

  gFull <- disimXdisimKernGradient(kern, kern, t, t2, covGrad)

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

  option <- 5
  l <- sqrt(2/disimKern1$inverseWidth)
  disimH1 <- disimComputeH(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l, option)
  h1 <- disimH1$h
  dh1_ddelta <- disimH1$dh_ddelta
  dh1_dD1 <- disimH1$dh_dD_j
  dh1_dD2 <- disimH1$dh_dD_k
  dh1_dl <- disimH1$dh_dl

  disimH2 <- disimComputeH(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l, option)
  h2 <- disimH2$h
  dh2_ddelta <- disimH2$dh_ddelta
  dh2_dD2 <- disimH2$dh_dD_j
  dh2_dD1 <- disimH2$dh_dD_k
  dh2_dl <- disimH2$dh_dl

  disimHp1 <- disimComputeHPrime(t1, t2, disimKern1$di_decay, disimKern1$decay, disimKern2$decay, l, option)

  hp1 <- disimHp1$h
  dhp1_ddelta <- disimHp1$dh_ddelta
  dhp1_dD1 <- disimHp1$dh_dD_j
  dhp1_dD2 <- disimHp1$dh_dD_k
  dhp1_dl <- disimHp1$dh_dl

  disimHp2 <- disimComputeHPrime(t2, t1, disimKern1$di_decay, disimKern2$decay, disimKern1$decay, l, option)
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
  K <- h1 + t(h2) + hp1 + t(hp2)
  K <- 0.5*K*sqrt(pi)
  var2 <- C0*C1*C2

  dk_ddelta <- sum(covGrad*dK_ddelta)*0.5*sqrt(pi)*l*var2
  dk_dD1 <- sum(covGrad*dK_dD1)*0.5*sqrt(pi)*l*var2
  dk_dD2 <- sum(covGrad*dK_dD2)*0.5*sqrt(pi)*l*var2
  dk_dl <- sum(covGrad*(dK_dl*0.5*sqrt(pi)*l + K))*var2
  K <- l*K
  dk_dC0 <- C1*C2*sum(covGrad*K)
  dk_dC1 <- C0*C2*sum(covGrad*K)
  dk_dC2 <- C0*C1*sum(covGrad*K)

  dk_dDIVariance <- dk_dC0
  dk_dDisim1Variance <- dk_dC1*0.5/C1
  dk_dDisim2Variance <- dk_dC2*0.5/C2

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern1$inverseWidth*sqrt(disimKern1$inverseWidth))*dk_dl

  K <- var2*K

  g1 <- c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD1, dk_dDisim1Variance)
  g2 <- c(0, 0, 0, dk_dD2, dk_dDisim2Variance)

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
  C_j <- sqrt(rbfKern$variance)
  D_i <- disimKern$decay
  delta <-  disimKern$di_decay

  invLDiffT <- 1/l*diffT
  halfLD_i <- 0.5*l*D_i
  halfLDelta <- 0.5*l*delta

  prefact <- C_0 * C_i * C_j * sqrt(pi)/2 * l

  lnCommon <- - complexLog(delta - D_i)
  lnPart1 <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)
  lnPart2 <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)

  lnFact1 <- halfLDelta^2 - delta * diffT
  lnFact2 <- halfLD_i^2 - D_i * diffT

  gradln1 <- gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, l/2, l/2)
  dlnPart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2 <-gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, l/2, l/2)
  dlnPart2 <- gradln2$dlnPart
  m2 <- gradln2$m

  dK_dD <- k * (1/(delta-D_i)) + prefact * ((l*halfLD_i - diffT) * exp(lnCommon + lnFact2 + lnPart2) + dlnPart2 * exp(lnCommon + lnFact2 - m2))
  dk_dD <- sum(dK_dD*covGrad)

  dK_ddelta <- k * (-1/(delta-D_i)) + prefact * ((l*halfLDelta - diffT) * exp(lnCommon + lnFact1 + lnPart1) + dlnPart1 * exp(lnCommon + lnFact1 - m1))
  dk_ddelta <- sum(dK_ddelta*covGrad)

  gradln1 <- gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, delta/2 + invLDiffT/l, delta/2 - t2Mat/l2)
  dlnPart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2 <-gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, D_i/2 - t2Mat/l2, D_i/2 + invLDiffT/l);
  dlnPart2 <- gradln2$dlnPart
  m2 <- gradln2$m

  dK_dl <- k/l + prefact * (delta*halfLDelta * exp(lnCommon + lnFact1 + lnPart1) + D_i*halfLD_i * exp(lnCommon + lnFact2 + lnPart2) + dlnPart1 * exp(lnCommon + lnFact1 - m1) + dlnPart2 * exp(lnCommon + lnFact2 - m2))
  dk_dl <- sum(dK_dl*covGrad)

  dk_dC_i <- sum(k*covGrad)/C_i
  dk_dC_0 <- sum(k*covGrad)/C_0
  dk_dRbfVariance <- 0.5*sum(k*covGrad)/rbfKern$variance

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern$inverseWidth*sqrt(disimKern$inverseWidth))*dk_dl
  dk_dDisimVariance <- dk_dC_i*0.5/C_i
  dk_dDIVariance <- dk_dC_0*0.5/C_0

  g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dDisimVariance))
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

  lnCommon1 <- - complexLog(2*delta) -delta * t2Mat - D_i * t1Mat + halfLDelta^2

  lnFact1 <- complexLog(2 * delta) - complexLog(delta^2 - D_i^2)
  lnPart1 <- lnDiffErfs(halfLDelta - t2Mat/l, halfLDelta)

  lnFact2 <- (D_i - delta) * t1Mat - complexLog(delta - D_i)
  lnPart2a <- lnDiffErfs(halfLDelta, halfLDelta - t1Mat/l)
  lnPart2b <- lnDiffErfs(halfLDelta, halfLDelta - t2Mat/l)

  lnFact3 <- (D_i + delta) * t1Mat - complexLog(delta + D_i)
  lnPart3 <- lnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT)

  lnFact4 <- 2*delta*t2Mat + (D_i - delta) * t1Mat - complexLog(delta - D_i)
  lnPart4 <- lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l)

  lnCommon2 <- - complexLog(delta^2 - D_i^2) - delta * t2Mat - D_i * t1Mat + halfLD_i^2
  lnPart5 <- lnDiffErfs(halfLD_i - t1Mat/l, halfLD_i)

  lnFact6 <- (D_i + delta) * t2Mat
  lnPart6 <- lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT)

  K1 <- exp( lnCommon1 + lnFact1 + lnPart1) +exp(lnCommon1 + lnFact2 + lnPart2a) +exp(lnCommon1 + lnFact2 + lnPart2b) +exp(lnCommon1 + lnFact3 + lnPart3) +exp(lnCommon1 + lnFact4 + lnPart4)

  K2 <- exp( lnCommon2 + lnPart5) +exp(lnCommon2 + lnFact6 + lnPart6)

  prefact <- 0.5*sqrt(pi)*l*disimKern$di_variance*sqrt(disimKern$variance)
  K <- prefact*Re(K1+K2)

  dcommon1 <- - 1/delta - t2Mat + l*halfLDelta
  dfact1 <- 1/delta - 2*delta/(delta^2 - D_i^2)
  dfact2 <- -t1Mat - 1/(delta - D_i)
  dfact3 <- t1Mat - 1/(delta + D_i)
  dfact4 <- 2*t2Mat - t1Mat - 1/(delta - D_i)
  dcommon2 <- - 2*delta/(delta^2 - D_i^2) - t2Mat
  dfact6 <- t2Mat

  gradln1 <- gradLnDiffErfs(halfLDelta - t2Mat/l, halfLDelta, l/2, l/2)
  dpart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2a <-gradLnDiffErfs(halfLDelta, halfLDelta - t1Mat/l, l/2, l/2)
  dpart2a <- gradln2a$dlnPart
  m2a <- gradln2a$m

  gradln2b <-gradLnDiffErfs(halfLDelta, halfLDelta - t2Mat/l, l/2, l/2)
  dpart2b <- gradln2b$dlnPart
  m2b <- gradln2b$m

  gradln3 <- gradLnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT, l/2, l/2)
  dpart3 <- gradln3$dlnPart
  m3 <- gradln3$m

  gradln4 <- gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, l/2, l/2)
  dpart4 <- gradln4$dlnPart
  m4 <- gradln4$m

  dK_ddelta <- prefact * (dcommon1 * K1 + dcommon2 * K2 + dfact1 * exp(lnCommon1 + lnFact1 + lnPart1) + dfact2 * exp(lnCommon1 + lnFact2 + lnPart2a) + dfact2 * exp(lnCommon1 + lnFact2 + lnPart2b) + dfact3 * exp(lnCommon1 + lnFact3 + lnPart3) + dfact4 * exp(lnCommon1 + lnFact4 + lnPart4) + dfact6 * exp(lnCommon2 + lnFact6 + lnPart6) + dpart1 * exp(lnCommon1 + lnFact1 - m1) + dpart2a * exp(lnCommon1 + lnFact2 - m2a) + dpart2b * exp(lnCommon1 + lnFact2 - m2b) + dpart3 * exp(lnCommon1 + lnFact3 - m3) + dpart4 * exp(lnCommon1 + lnFact4 - m4))

  dcommon1 <- delta * halfLDelta
  dcommon2 <- D_i * halfLD_i

  gradln1 <- gradLnDiffErfs(halfLDelta - t2Mat/l, halfLDelta, delta/2 + t2Mat/l2, delta/2)
  dpart1 <- gradln1$dlnPart
  m1 <- gradln1$m

  gradln2a <-gradLnDiffErfs(halfLDelta, halfLDelta - t1Mat/l, delta/2, delta/2 + t1Mat/l2)
  dpart2a <- gradln2a$dlnPart
  m2a <- gradln2a$m

  gradln2b <-gradLnDiffErfs(halfLDelta, halfLDelta - t2Mat/l, delta/2, delta/2 + t2Mat/l2)
  dpart2b <- gradln2b$dlnPart
  m2b <- gradln2b$m

  gradln3 <- gradLnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT, delta/2 - t1Mat/l2, delta/2 - invLDiffT/l)
  dpart3 <- gradln3$dlnPart
  m3 <- gradln3$m

  gradln4 <- gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, delta/2 + invLDiffT/l, delta/2 - t2Mat/l2)
  dpart4 <- gradln4$dlnPart
  m4 <- gradln4$m

  gradln5 <- gradLnDiffErfs(halfLD_i - t1Mat/l, halfLD_i, D_i/2 + t1Mat/l2, D_i/2)
  dpart5 <- gradln5$dlnPart
  m5 <- gradln5$m

  gradln6 <- gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, D_i/2 - t2Mat/l2, D_i/2 + invLDiffT/l)
  dpart6 <- gradln6$dlnPart
  m6 <- gradln6$m

  dK_dl <- Re(prefact * (dcommon1 * K1 + dcommon2 * K2 +dpart1 * exp(lnCommon1 + lnFact1 - m1) +dpart2a * exp(lnCommon1 + lnFact2 - m2a) +dpart2b * exp(lnCommon1 + lnFact2 - m2b) +dpart3 * exp(lnCommon1 + lnFact3 - m3) +dpart4 * exp(lnCommon1 + lnFact4 - m4) +dpart5 * exp(lnCommon2 - m5) +dpart6 * exp(lnCommon2 + lnFact6 - m6))) + K / l


  dcommon1 <- - t1Mat
  dfact1 <- 2 * D_i / (delta^2 - D_i^2)
  dfact2 <- t1Mat + 1/(delta - D_i)
  dfact3 <- t1Mat - 1/(delta + D_i)
  dfact4 <- t1Mat + 1/(delta - D_i)
  dcommon2 <- 2 * D_i / (delta^2 - D_i^2) - t1Mat + l*halfLD_i
  dfact6 <- t2Mat

  gradln5 <- gradLnDiffErfs(halfLD_i - t1Mat/l, halfLD_i, l/2, l/2)
  dpart5 <- gradln5$dlnPart
  m5 <- gradln5$m

  gradln6 <- gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, l/2, l/2)
  dpart6 <- gradln6$dlnPart
  m6 <- gradln6$m

  dK_dD = prefact * (dcommon1 * K1 + dcommon2 * K2 + dfact1 * exp(lnCommon1 + lnFact1 + lnPart1) + dfact2 * exp(lnCommon1 + lnFact2 + lnPart2a) + dfact2 * exp(lnCommon1 + lnFact2 + lnPart2b) + dfact3 * exp(lnCommon1 + lnFact3 + lnPart3) + dfact4 * exp(lnCommon1 + lnFact4 + lnPart4) + dfact6 * exp(lnCommon2 + lnFact6 + lnPart6) + dpart5 * exp(lnCommon2 - m5) + dpart6 * exp(lnCommon2 + lnFact6 - m6))

  dk_ddelta <- sum(dK_ddelta*covGrad)
  dk_dl <- sum(dK_dl*covGrad)
  dk_dD <- sum(dK_dD*covGrad)

  dk_dDIVariance <- sum(K*covGrad)/disimKern$di_variance
  dk_dSimVariance <- .5 * sum(K*covGrad)/disimKern$variance

  dk_dinvWidth <- -0.5*sqrt(2)/(disimKern$inverseWidth* sqrt(disimKern$inverseWidth))*dk_dl

  g1 <- Re(c(dk_ddelta, dk_dinvWidth, dk_dDIVariance, dk_dD, dk_dSimVariance))
  g2 <- c(0, 0, 0)

  g <- list(g1=g1, g2=g2)
  return (g)
}








