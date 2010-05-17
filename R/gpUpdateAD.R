gpUpdateAD <- function (model, X) {
  
  if (nargs() < 2)
    X = model$X

  if (model.approx == "ftc") {
    ## Compute the inner product values.
    if (!"S" %in% names(model)) {
      if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
	for (i in 1:model$d)
	  model$innerProducts[1, i] = model$m[, i] %*% model$invK_uu %*% model$m[, i]
      } else {
	for (i in 1:model$d) {
	  ind = gpDataIndices(model, i)
	  model$innerProducts[1, i] = model$m[, i] %*% model$invK_uu[[i]] %*% model$m[, i]
	}
      }
    }
  } else if (model$approx %in% c("dtc", "dtcvar")) {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      ## Compute A = invBetaK_uu + K_uf*K_uf'
      K_uf2 = model$K_uf %*% model$K_uf
      model$A = (1/model$beta)%*%model$K_uu + K_uf2
      ## This can become unstable when K_uf2 is low rank.
      invA = .jitCholInv(model.A, silent=TRUE)
      model$Ainv = invA$invM
      model$logDetA = 2* sum( log ( diag(invA$chol) ) )

      ## compute inner products
      for (i in 1:model$d) {
	E = model$K_uf %*% model$m[, i]
	model$innerProducts[1, i] = model$beta %*%
	  (model$m[:, i]%*%model$m[, i] - E%*%model$Ainv%*%E)
      }

      if (model$approx == "dtcvar") {
	model$diagD = model$beta %*% (model$diagK
	  - colSums(model$K_uf * (model$invK_uu%*%model$K_uf)))
      }
    } else {
      if (!model$isMissingData)
	K_uf2 = model$K_uf%*%model$K_uf

      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	## Compute A = invBetaK_uu + K_uf*K_uf'
	if (model$isMissingData)
	  K_uf2 = model$K_uf[, ind]%*%model$K_uf[, ind]

	model$A[[i]] = (1/model$beta)%*%model$K_uu + K_uf2
	## This can become unstable when K_uf2 is low rank.
	invA = .jitCholInv(model.A[[i]], silent=TRUE)
	model$Ainv[[i]] = invA$invM
	model$logDetA[i] = 2* sum( log ( diag(invA$chol) ) )
	## compute inner products
	E = model$K_uf[, ind]%*%model$m[ind, i]
	model$innerProducts[1, i] = 
	  model$beta%*%(model$m[ind, i]%*%model$m[ind, i] - E%*%model$Ainv[[i]]%*%E)
      }

      if (model$approx == "dtcvar")
	stop("Non spherical implementation for dtcvar not yet done.")
    }
  } else if (model.approx == "fitc") {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      model$A = 1/model$beta%*%model$K_uu
      K_ufDinvm = matrix(0,model$k, model$d)
      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	model$D[[i]] = diag(length(ind)) + model$beta%*%model$K[[i]] -
	  model$beta%*%model$K_uf[, ind]%*%model$invK_uu%*%model$K_uf[, ind]
	invD = .jitCholInv(model.D[[i]], silent=TRUE)
	model$Dinv[[i]] = invD$invM
	model$logDetD[i] = 2* sum( log ( diag(invD$chol) ) )
	K_ufDinvK_uf = model.$K_uf[, ind]%*%model$Dinv[[i]]%*%model$K_uf[, ind]
	model$A = model$A + K_ufDinvK_uf
	Dinvm[[i]] = model$Dinv{i}%*%model$m[ind, ]
	K_ufDinvm = K_ufDinvm + model$K_uf[, ind]%*%Dinvm[[i]]
      }
      ## This can become unstable when K_ufDinvK_uf is low rank.
      invA = .jitCholInv(model.A, silent=TRUE)
      model$Ainv = invA$invM
      model$logDetA = 2* sum( log ( diag(invA$chol) ) )
      ## compute inner products
      for (i in 1:model$d)
	model$innerProducts[1, i] = - model$beta%*%K_ufDinvm[, i]%*%model$Ainv
	  %*%K_ufDinvm[, i]

      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	for (j in 1:model$d)
	  model$innerProducts[1, j] = model$innerProducts[1, j]
	    + model$beta%*%Dinvm[[i]][, j]%*%model$m[ind, j]
      }
    } else {
      for (j in 1:model$d) {
	model$A[[j]] = 1/model$beta%*%model$K_uu
	K_ufDinvm = matrix(0, model$k, model$d)
	for (i in 1:length(model$blockEnd)) {
	  ind = gpDataIndices(model, j, i)
	  model$D[[i, j]] = diag(length(ind)) + model$beta%*%model$K[[i, j]] -
	      model$beta%*%model$K_uf[, ind]%*%model$invK_uu%*%model$K_uf[, ind]
	  invD = .jitCholInv(model.D[[i,j]], silent=TRUE)
	  model$Dinv[[i,j]] = invD$invM
	  model$logDetD[i,j] = 2* sum( log ( diag(invD$chol) ) )
	  K_ufDinvK_uf = model$K_uf[, ind]%*%model$Dinv[[i, j]]%*%model$K_uf[, ind]
	  model$A[[j]] = model$A[[j]] + K_ufDinvK_uf
	  Dinvm[[i]][ind, j] = model$Dinv[[i, j]]%*%model$m[ind, j]
	  K_ufDinvm[, j] = K_ufDinvm[, j] + model$K_uf[, ind]%*%Dinvm[[i]][ind, j]
	}
	## This can become unstable when K_ufDinvK_uf is low rank.
	invA = .jitCholInv(model.A[[j]], silent=TRUE)
	model$Ainv[[i]] = invA$invM
	model$logDetA[j] = 2* sum( log ( diag(invA$chol) ) )
      }
      
      ## compute inner products
      for (j = 1:model$d)
	model$innerProducts[1, j] = - model$beta%*%K_ufDinvm[, j]%*%model$Ainv[[j]]
	  %*%K_ufDinvm[, j]

      for (i in 1:length(model$blockEnd)) {
	for (j in 1:model$d) {
	  ind = gpDataIndices(model, j, i)
	  model$innerProducts[1, j] = model$innerProducts[1, j]
	    + model$beta%*%Dinvm[[i]][ind, j]%*%model$m[ind, j]
	}
      }
    }
  } else
    stop("Unknown approximating criterion.")

  return (model)
}