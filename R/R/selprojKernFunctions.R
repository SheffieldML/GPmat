#
# selproj kernel
#   wrap around another simple kernel
#   returning its value if x[1,]==kern$options$expmask,
#   or x[1,]==0 or kern$options$expmask==0
#


selprojKernDiagCompute <- function (kern, x) {
  d <- dim(as.array(x))[1]
  k <- array(0, dim=d)

  if (dim(as.array(x))[2] > 1) {
    if (kern$expmask > 0)
      I <- (x[,1] == kern$expmask)
    else
      I <- 1:d

    k[I] <- kernDiagCompute(kern$comp[[1]], x[I,-1])
  }

  return (k)
}



selprojKernExtractParam <- function (kern, only.values=TRUE,
                                     untransformed.values=FALSE) {
  funcName <- paste(kern$comp[[1]]$type, "KernExtractParam", sep="")
  func <- get(funcName, mode="function")

  params <- func(kern$comp[[1]], only.values=only.values,
                 untransformed.values=untransformed.values)

  if (!only.values) {
    names(params) <- paste(kern$comp[[1]]$type, names(params), sep='_')
  }

  return (params)
}



selprojKernExpandParam <- function (kern, params) {
  funcName <- paste(kern$comp[[1]]$type, "KernExpandParam", sep="")
  func <- get(funcName, mode="function")
  kern$comp[[1]] <- func(kern$comp[[1]], params)
  
  return (kern)
}



selprojKernCompute <- function (kern, x, x2) {
  if ( nargs()>2 ) {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]
    k <- array(0, dim=c(d1, d2))

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1) {
      if (kern$expmask > 0) {
        I1 <- (x[,1] == kern$expmask)
        I2 <- (x2[,1] == kern$expmask)

        if (!any(I1) || !any(I2))
          return (k)
      } else {
        I1 <- 1:d1
        I2 <- 1:d2
      }

      k[I1,I2] <- kernCompute(kern$comp[[1]], x[I1,-1], x2[I2,-1])
    }
  } else {
    d1 <- dim(as.array(x))[1]
    k <- array(0, dim=c(d1, d1))

    if (dim(as.array(x))[2] > 1) {
      if (kern$expmask > 0) {
        I1 <- (x[,1] == kern$expmask)

        if (!any(I1))
          return (k)
      } else {
        I1 <- 1:d1
      }
      
      k[I1,I1] <- kernCompute(kern$comp[[1]], x[I1,-1])
    }
  }
  return (k)  
}



selprojXselprojKernCompute <- function (kern1, kern2, x, x2) {
  funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernCompute", sep="")
  func <- get(funcName, mode="function")
  
  if ( nargs()>3 ) {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]
    k <- array(0, dim=c(d1, d2))

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1) {
      # Zero output if different non-zero expmasks
      if (kern1$expmask <= 0 || kern2$expmask <= 0 ||
          kern1$expmask == kern2$expmask) {
        expmask <- max(kern1$expmask, kern2$expmask)
        if (expmask > 0) {
          I1 <- (x[,1] == expmask)
          I2 <- (x2[,1] == expmask)

          if (!any(I1) || !any(I2))
            return (k)
        } else {
          I1 <- 1:d1
          I2 <- 1:d2
        }

        k[I1,I2] <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-1], x2[I2,-1])
      }
    }
  } else {
    d1 <- dim(as.array(x))[1]
    k <- array(0, dim=c(d1, d1))

    if (dim(as.array(x))[2] > 1) {
      if (kern1$expmask <= 0 || kern2$expmask <= 0 ||
          kern1$expmask == kern2$expmask) {
        expmask <- max(kern1$expmask, kern2$expmask)
        if (expmask > 0) {
          I1 <- (x[,1] == expmask)
          
        if (!any(I1))
          return (k)
        } else {
          I1 <- 1:d1
        }
      
        k[I1,I1] <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-1])
      }
    }
  }
  return (k)  
}



selprojKernParamInit <- function (kern) {
  kern <- cmpndKernParamInit(kern)

  for (i in seq_along(kern$comp))
    if (kern$comp[[i]]$type == "multi")
      error("selproj kernel does not support wrapping a multi kernel.")
  
  kern$expmask <- kern$options$expmask
  kern$transforms <- kern$comp[[1]]$transforms
  if ("transformArgs" %in% names(kern$comp[[1]]))
    kern$transformArgs <- kern$comp[[1]]$transformArgs
  
  return (kern)
}



selprojKernGradient <- function (kern, x, x2, covGrad) {
  funcName <- paste(kern$comp[[1]]$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")

  if ( nargs()<4 ) {
    covGrad <- x2

    d1 <- dim(as.array(x))[1]

    if (dim(as.array(x))[2] > 1) {
      if (kern$expmask > 0) {
        I1 <- (x[,1] == kern$expmask)

        if (!any(I1))
          return (rep(0, kern$nParams))
      } else {
        I1 <- 1:d1
      }
    }

    g <- func(kern$comp[[1]], x[I1,-1], covGrad[I1,I1])
  } else {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1) {
      if (kern$expmask > 0) {
        I1 <- (x[,1] == kern$expmask)
        I2 <- (x2[,1] == kern$expmask)

        if (!any(I1) || !any(I2))
          return (rep(0, kern$nParams))
      } else {
        I1 <- 1:d1
        I2 <- 1:d2
      }
    }
    
    g <- func(kern$comp[[1]], x[I1,-1], x2[I2,-1], covGrad[I1,I2])
  }

  return (g)
}


selprojXselprojKernGradient <- function (kern1, kern2, x, x2, covGrad) {
  funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")

  if ( nargs()<5 ) {
    covGrad <- x2

    d1 <- dim(as.array(x))[1]

    if (dim(as.array(x))[2] > 1 &&
        (kern1$expmask <= 0 || kern2$expmask <= 0 ||
         kern1$expmask == kern2$expmask)) {
      expmask <- max(kern1$expmask, kern2$expmask)

      if (expmask > 0) {
        I1 <- (x[,1] == expmask)

        return (list(g1=rep(0, kern1$nParams), g2=rep(0, kern2$nParams)))
      } else {
        I1 <- 1:d1
      }

      g <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-1], covGrad[I1,I1])
    } else {
      g <- list(g1=rep(0, kern1$nParams), g2=rep(0, kern2$nParams))
    }
  } else {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1 &&
        (kern1$expmask <= 0 || kern2$expmask <= 0 ||
         kern1$expmask == kern2$expmask)) {
      expmask <- max(kern1$expmask, kern2$expmask)

      if (expmask > 0) {
        I1 <- (x[,1] == expmask)
        I2 <- (x2[,1] == expmask)

        if (!any(I1) || !any(I2))
          return (list(g1=rep(0, kern1$nParams), g2=rep(0, kern2$nParams)))
      } else {
        I1 <- 1:d1
        I2 <- 1:d2
      }

      g <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-1], x2[I2,-1], covGrad[I1,I2])
    } else {
      g <- list(g1=rep(0, kern1$nParams), g2=rep(0, kern2$nParams))
    }
  }

  return (g)
}


selprojKernDisplay <- function (kern, spaceNum=0) {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Selective projection kernel:\n")
  cat(spacing)
  cat("Experiment mask:", kern$expmask, "\n")
  for(i in seq(along=kern$comp)) 
    kernDisplay(kern$comp[[i]], spaceNum+2)
}
