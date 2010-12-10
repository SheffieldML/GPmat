gpUpdateKernels <- function (model, X, X_u) {
  jitter =  1e-6

  if (model$approx == "ftc") {
    ## Long term should allow different kernels in each dimension here.
    model$K_uu = kernCompute(model$kern, X)

    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      ## Add inverse beta to diagonal if it exists.
      if ("beta" %in% names(model) && length(model$beta)>0) {
	model$K_uu[seq(1,length(model$K_uu),by= dim(model$K_uu)[1]+1)] = 
	  model$K_uu[seq(1,length(model$K_uu),by= dim(model$K_uu)[1]+1)] + 1/model$beta
      }
      invK = .jitCholInv(model$K_uu, silent=TRUE) ## pdinv + jitChol combined
      model$invK_uu = invK$invM
      model$logDetK_uu = 2* sum( log ( diag(invK$chol) ) )
    } else {
      model$invK_uu=list(); model$logDetK_uu=matrix(0,1,model$d)
      for (i in 1:model$d) {
	if ("beta" %in% names(model) && length(model$beta)>0) {
	  if (dim(as.matrix(model$beta))[2] == model$d)
	    betaAdd = model$beta[, i]
	  else
	    betaAdd = model$beta

	  model$K_uu[seq(1,length(model$K_uu),by= dim(model$K_uu)[1]+1)] = 
	    model$K_uu[seq(1,length(model$K_uu),by= dim(model$K_uu)[1]+1)] + 1/betaAdd
	}
	ind = gpDataIndices(model, i)
	invK = .jitCholInv(model$K_uu[ind,ind], silent=TRUE) ## pdinv + jitChol combined
	model$invK_uu[[i]] = invK$invM
	model$logDetK_uu[i] = 2* sum( log ( diag(invK$chol) ) )
      }
    }
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    model$K_uu = kernCompute(model$kern, X_u)

    if ((!"whiteVariance" %in% names(model$kern)) || model$kern$whiteVariance == 0) {
      ## There is no white noise term so add some jitter.
      model$K_uu = model$K_uu + diag.spam(jitter, dim(model$K_uu)[1]) ## need 'spam'
      #sparseDiag(matrix(jitter, dim(model$K_uu)[1], 1))
    }
    model$K_uf = kernCompute(model$kern, X_u, X)
    invK = .jitCholInv(model$K_uu, silent=TRUE) ## pdinv + jitChol combined
    model$invK_uu = invK$invM
    model$logDetK_uu = 2* sum( log ( diag(invK$chol) ) )
  }

  if (model$approx %in% c("dtcvar", "fitc"))
    model$diagK = kernDiagCompute(model$kern, X)
  else if (model$approx == "pitc") {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      model$K=list()
      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	model$K[[i]] = kernCompute(model$kern, X[ind, ,drop=FALSE])
      }
    } else {
      model$K = matrix(0, length(model$blockEnd), model$d)
      model$K = lapply(split(model$K,row(model$K)), split, 1:model$d)
      for (j in 1:model$d) {
	for (i in 1:length(model$blockEnd)) {
	  ind = gpDataIndices(model, j, i)
	  model$K[[i]][[j]] = kernCompute(model$kern, X[ind, ,drop=FALSE])
	}
      }
    }
  }

  model = gpUpdateAD(model, X)

  return (model)
}