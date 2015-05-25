gpUpdateAD <- function (model, X=model$X) {

  model$beta = drop(model$beta)
  if (model$approx == "ftc") {
    ## Compute the inner product values.
    if (!"S" %in% names(model)) {
      model$innerProducts = matrix(0, 1, model$d)
      if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
	for (i in 1:model$d) {
	  model$innerProducts[1, i] = t(model$m[, i,drop=FALSE])%*%model$invK_uu%*%model$m[, i,drop=FALSE]
	}
      } else {
	for (i in 1:model$d) {
	  ind = gpDataIndices(model, i)
	  model$innerProducts[1, i] = t(model$m[ind, i,drop=FALSE])%*%model$invK_uu[[i]]%*%model$m[ind, i,drop=FALSE]
	}
      }
    }
  } else if (model$approx %in% c("dtc", "dtcvar")) {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      ## Compute A = invBetaK_uu + K_uf%*%K_uf'
      K_uf2 = model$K_uf %*% t(model$K_uf)
      model$A = 1/model$beta * model$K_uu + K_uf2
      ## This can become unstable when K_uf2 is low rank.
      invA = .jitCholInv(model$A, silent=TRUE)
      model$Ainv = invA$invM
      model$logDetA = 2* sum(log(diag(invA$chol)))

      ## compute inner products
      model$innerProducts = matrix(0, 1, model$d)
      for (i in 1:model$d) {
	E = model$K_uf %*% model$m[, i]
	model$innerProducts[1, i] = model$beta *
	  (t(model$m[, i,drop=FALSE])%*%model$m[, i,drop=FALSE] - t(E)%*%model$Ainv%*%E)
      }

      if (model$approx == "dtcvar") {
	model$diagD = model$beta *
	  (model$diagK - colSums(model$K_uf * (model$invK_uu%*%model$K_uf)))
      }
    } else {
      model$A=list(); model$logDetA=matrix(0,1,model$d)
      if (!model$isMissingData)
	K_uf2 = model$K_uf%*%t(model$K_uf)

      model$innerProducts = matrix(0, 1, model$d)
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	## Compute A = invBetaK_uu + K_uf%*%K_uf'
	if (model$isMissingData)
	  K_uf2 = model$K_uf[, ind,drop=FALSE]%*%t(model$K_uf[, ind,drop=FALSE])

	model$A[[i]]= (1/model$beta) * model$K_uu+K_uf2
	## This can become unstable when K_uf2 is low rank.
	invA = .jitCholInv(model$A[[i]], silent=TRUE)
	model$Ainv[[i]] = invA$invM
	model$logDetA[i] = 2* sum( log ( diag(invA$chol) ) )
	## compute inner products
	E = model$K_uf[, ind,drop=FALSE]%*%model$m[ind, i,drop=FALSE]
	model$innerProducts[1, i] = (model$beta) * 
	  (t(model$m[ind, i,drop=FALSE])%*%model$m[ind, i,drop=FALSE] - t(E)%*%model$Ainv[[i]]%*%E)
      }

      if (model$approx == "dtcvar")
	stop("Non spherical implementation for dtcvar not yet done.")
    }
  } else if (model$approx == "fitc") {
    model$L = t(.jitChol(model$K_uu)$chol)

    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      model$diagD = 1 + (model$beta)*model$diagK
	  - model$beta*t(colSums(model$K_uf * (model$invK_uu%*%model$K_uf)))

      model$Dinv = diag.spam(drop(1/model$diagD)) #sparseDiag(1/model$diagD)
      K_ufDinvK_uf = model$K_uf%*%model$Dinv%*%t(model$K_uf)
      model$A = 1/model$beta*model$K_uu + K_ufDinvK_uf
      ## This can become unstable when K_ufDinvK_uf is low rank.
      invA = .jitCholInv(model$A, silent=TRUE)
      model$Ainv = invA$invM
      model$logDetA = 2*sum(log(diag(invA$chol)))
      model$detDiff = - log(model$beta)*model$k +
	log(det(diag(model$k) + model$beta*K_ufDinvK_uf%*%model$invK_uu))
      ## compute inner products
      model$innerProducts = matrix(0, 1, model$d)
      for (i in 1:model$d) {
	Dinvm = model$Dinv %*% model$m[, i]
	K_ufDinvm = model$K_uf%*%Dinvm
	model$innerProducts[1, i] = model$beta *
	    (t(Dinvm)%*%model$m[, i,drop=FALSE] - t(K_ufDinvm)%*%model$Ainv%*%K_ufDinvm)
      }
    
      ## Computations from Ed's implementation.
      model$V = solve(model$L, model$K_uf) #model$L \ model$K_uf
      model$V = model$V / kronecker(matrix(1,model$k,1), t(sqrt(model$diagD))) #repmat(sqrt(model$diagD)', model$k, 1)
      model$Am = 1/model$beta*diag(model$k) + model$V%*%t(model$V)
      model$Lm = t(.jitChol(model$Am)$chol)
      model$invLmV = solve(model$Lm, model$V)
      model$scaledM = model$m / kronecker(matrix(1,1,model$d), sqrt(model$diagD))
      model$bet = model$invLmV%*%model$scaledM
    } else {
      model$innerProducts = matrix(0, 1, model$d)
      model$logDetA=matrix(0,1,model$d); model$detDiff=matrix(0,1,model$d)
      model$diagD=list(); model$Dinv=list(); model$A=list(); model$Ainv=list()
      model$V=list(); model$Am=list(); model$Lm=list(); model$invLmV=list()
      model$scaledM=list(); model$bet=list()
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$diagD[[i]] = 1 + model$beta*model$diagK[ind]
	  - model$beta*t(colSums(model$K_uf[, ind,drop=FALSE] * (model$invK_uu%*%model$K_uf[, ind,drop=FALSE])))
	model$Dinv[[i]] = diag.spam(drop(1/model$diagD[[i]]))
	K_ufDinvK_uf = model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]%*%t(model$K_uf[, ind,drop=FALSE])
	model$A[[i]] = 1 / model$beta*model$K_uu + K_ufDinvK_uf
	## This can become unstable when K_ufDinvK_uf is low rank.
	invA = .jitCholInv(model$A[[i]], silent=TRUE)
	model$Ainv[[i]] = invA$invM
	model$logDetA[i] = 2*sum(log(diag(invA$chol)))
	model$detDiff[i] = - log(model$beta)*model$k
	  + log(det(diag(model$k) + model$beta*K_ufDinvK_uf%*%model$invK_uu))

	## compute inner products
	Dinvm = model$Dinv[[i]]%*%model$m[ind, i,drop=FALSE]
	K_ufDinvm = model$K_uf[, ind,drop=FALSE]%*%Dinvm
	model$innerProducts[1, i] = model$beta*(t(Dinvm)%*%model$m[ind, i,drop=FALSE] - t(K_ufDinvm)%*%model$Ainv[[i]]%*%K_ufDinvm)

	## Computations from Ed's implementation.
	model$V[[i]] = solve(model$L, model$K_uf[, ind,drop=FALSE])
	model$V[[i]] = model$V[[i]] / kronecker(matrix(1,model$k,1), t(sqrt(model$diagD[[i]])))
	model$Am[[i]] = 1/model$beta * diag(model$k) + model$V[[i]]%*%t(model$V[[i]])
	model$Lm[[i]] = t(.jitChol(model$Am[[i]])$chol)
	model$invLmV[[i]] = solve(model$Lm[[i]], model$V[[i]])
	model$scaledM[[i]] = model$m[ind, i,drop=FALSE] / sqrt(model$diagD[[i]])
	model$bet[[i]] = model$invLmV[[i]]%*%model$scaledM[[i]]
      }
    }
  } else if (model$approx == "pitc") {
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      model$A = 1/model$beta*model$K_uu
      K_ufDinvm = matrix(0,model$k, model$d)
      model$logDetD = matrix(0, 1, length(model$blockEnd))
      Dinvm = list(); model$Dinv = list(); model$D = list()
      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	model$D[[i]] = diag(length(ind)) + model$beta*model$K[[i]] -
	  model$beta*t(model$K_uf[, ind,drop=FALSE])%*%model$invK_uu%*%model$K_uf[, ind,drop=FALSE]
	invD = .jitCholInv(model$D[[i]], silent=TRUE)
	model$Dinv[[i]] = invD$invM
	model$logDetD[i] = 2* sum( log ( diag(invD$chol) ) )
	K_ufDinvK_uf = model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]%*%t(model$K_uf[, ind,drop=FALSE])
	model$A = model$A + K_ufDinvK_uf
	Dinvm[[i]] = model$Dinv[[i]]%*%model$m[ind, ,drop=FALSE]
	K_ufDinvm = K_ufDinvm + model$K_uf[, ind,drop=FALSE]%*%Dinvm[[i]]
      }
      ## This can become unstable when K_ufDinvK_uf is low rank.
      invA = .jitCholInv(model$A, silent=TRUE)
      model$Ainv = invA$invM
      model$logDetA = 2* sum(log(diag(invA$chol)))
      ## compute inner products
      model$innerProducts = matrix(0, 1, model$d)
      for (i in 1:model$d)
	model$innerProducts[1, i] = - model$beta*t(K_ufDinvm[, i,drop=FALSE])%*%model$Ainv%*%K_ufDinvm[, i,drop=FALSE]

      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	for (j in 1:model$d)
	  model$innerProducts[1, j] = model$innerProducts[1, j]
	    + model$beta*t(Dinvm[[i]][, j,drop=FALSE])%*%model$m[ind, j,drop=FALSE]
      }
    } else {
      model$A = list(); model$Ainv = list()
      model$D = matrix(0, length(model$blockEnd), model$d)
      model$D = lapply(split(model$D,row(model$D)), split, 1:model$d)
      model$Dinv = model$D
      model$logDetD = matrix(0, length(model$blockEnd), model$d)
      model$logDetA = matrix(0, 1, model$d)
      Dinvm = as.list(matrix(0,1,length(model$blockEnd)))
      Dinvm = lapply(Dinvm, function(x) x=matrix(0,model$N,model$d))
      for (j in 1:model$d) {
	model$A[[j]] = 1/model$beta*model$K_uu
	K_ufDinvm = matrix(0, model$k, model$d)
	for (i in 1:length(model$blockEnd)) {
	  ind = gpDataIndices(model, j, i)
	  model$D[[i]][[j]] = diag(length(ind)) + model$beta*model$K[[i]][[j]] -
	      model$beta*t(model$K_uf[, ind,drop=FALSE])%*%model$invK_uu%*%model$K_uf[, ind,drop=FALSE]
	  invD = .jitCholInv(model$D[[i]][[j]], silent=TRUE)
	  model$Dinv[[i]][[j]] = invD$invM
	  model$logDetD[i,j] = 2* sum( log ( diag(invD$chol) ) )
	  K_ufDinvK_uf = model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]][[j]]%*%t(model$K_uf[, ind,drop=FALSE])
	  model$A[[j]] = model$A[[j]] + K_ufDinvK_uf
	  Dinvm[[i]][ind, j] = model$Dinv[[i]][[j]]%*%model$m[ind, j,drop=FALSE]
	  K_ufDinvm[, j] = K_ufDinvm[, j,drop=FALSE] + model$K_uf[, ind,drop=FALSE]%*%Dinvm[[i]][ind, j,drop=FALSE]
	}
	## This can become unstable when K_ufDinvK_uf is low rank.
	invA = .jitCholInv(model$A[[j]], silent=TRUE)
	model$Ainv[[j]] = invA$invM
	model$logDetA[j] = 2* sum( log ( diag(invA$chol) ) )
      }
      
      model$innerProducts = matrix(0, 1, model$d)
      ## compute inner products
      for (j in 1:model$d)
	model$innerProducts[1, j] = - model$beta*t(K_ufDinvm[, j,drop=FALSE])%*%model$Ainv[[j]]%*%K_ufDinvm[, j,drop=FALSE]

      for (i in 1:length(model$blockEnd)) {
	for (j in 1:model$d) {
	  ind = gpDataIndices(model, j, i)
	  model$innerProducts[1, j] = model$innerProducts[1, j]
				    + model$beta*t(Dinvm[[i]][ind, j,drop=FALSE])%*%model$m[ind, j,drop=FALSE]
	}
      }
    }
  } else
    stop("Unknown approximating criterion.")

  return (model)
}
