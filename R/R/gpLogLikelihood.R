gpLogLikelihood <- function(model) {

  if (model$approx == "ftc") {
    ## No approximation, just do a full computation on K.
    ## For very high D, we use the matrix S which is M%*%M'
    if ("S" %in% names(model)) {
      ll = -0.5*(model$d*model$logDetK_uu + sum(model$invK_uu * model$S))
      return (ll)
    }
    ll = 0
    for (i in 1:dim(model$m)[2]) {
      if ((!"isSpherical" %in% names(model)) || model$isSpherical)
	ll = ll -.5*model$logDetK_uu - .5*t(model$m[, i,drop=FALSE])%*%model$invK_uu%*%model$m[, i,drop=FALSE]
      else {
	if (model$isMissingData)
	  m = model$m[model$indexPresent[[i]], i]
	else
	  m = model$m[, i]

	ll = ll - .5*model$logDetK_uu[i] - .5*t(m)%*%model$invK_uu[[i]]%*%m
      }
    }
  } else if (model$approx %in% c("dtc", "dtcvar")) {
    ## Deterministic training conditional
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      E = model$K_uf%*%model$m
      EET = E %*% t(E)
      if (length(model$beta)==1) {
	ll =  -0.5*(model$d*(-(model$N-model$k)*log(model$beta)
	      - model$logDetK_uu + model$logDetA) - (sum(model$Ainv*EET)
	      -sum(model$m * model$m))*model$beta)
	if (model$approx == "dtcvar")
	  ll = ll - model$d * 0.5*sum(model$diagD)
      } else
	stop("Not implemented variable length beta yet.")
    } else {
      ll = 0
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	e = model$K_uf[, ind,drop=FALSE]%*%model$m[ind, i,drop=FALSE]
	if (length(model$beta)==1) {
	  ll = ll - 0.5*((-(model$N-model$k)*log(model$beta)
		  - model$logDetK_uu + model$logDetA[i]) - (t(e)%*%model$Ainv[[i]]%*%e
		  - t(model$m[ind, i,drop=FALSE])%*%model$m[ind, i,drop=FALSE])*model$beta)
	  if(is.nan(ll))
	    stop("Log likelihood is NaN")

	  if (model$approx == "dtcvar")
	    stop("Not implemented dtcvar for non-spherical yet.")

	} else
	  stop("Not implemented variable length beta yet.")
      }
    }
  } else if (model$approx == "fitc") {
    ## Fully independent training conditional.
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      if (length(model$beta)==1) {
	if (FALSE) { ##  ಠ_ಠ ?
	  ## This is the original objective
	  Dinvm = model$Dinv %*% model$m
	  K_ufDinvm = model$K_uf %*% Dinvm
	  ll = -0.5*(model$d * (sum(log(model$diagD))
		-(model$N-model$k)*log(model$beta) + model$detDiff)
		+ (sum(Dinvm * model$m)
		- sum((model$Ainv%*%K_ufDinvm) * K_ufDinvm))*model$beta)
	  ll = ll - 0.5*model$N*model$d*log(2*pi)
	} else {
	  ## This is objective to match Ed Snelson's code
	  ll =  - model$d*(sum(log(diag(model$Lm)))
		+ 0.5*(-(model$N - model$k)*log(model$beta)
		+(model$N*log(2*pi) + sum(log(model$diagD)))))
	  for (i in 1:model$d)
	    ll = ll - 0.5*model$beta*(t(model$scaledM[, i,drop=FALSE])%*%model$scaledM[, i,drop=FALSE]
		    - t(model$bet[, i,drop=FALSE])%*%model$bet[, i,drop=FALSE])
	}
      } else
	stop("Variable length Beta not implemented yet.")
    } else {
      if (length(model$beta)==1) {
	if (FALSE) {
	  ll = 0
	  for (i in 1:model$d) {
	    ind = gpDataIndices(model, i)
	    Dinvm = model$Dinv[[i]]%*%model$m[ind, i,drop=FALSE]
	    K_ufDinvm = model$K_uf[, ind,drop=FALSE]%*%Dinvm
	    ll = ll -0.5*(sum(log(model$diagD[[i]]))
		    - (length(ind) - model$k)*log(model$beta)
		    + model$detDiff[i] + (sum(Dinvm * model$m[ind, i,drop=FALSE])
		    - sum((model$Ainv[[i]]%*%K_ufDinvm) * K_ufDinvm))*model$beta
		    + length(ind)*log(2*pi))
	  }
	} else {
	  ## This is objective to match Ed Snelson's code
	  ll = 0
	  for (i in 1:model$d) {
	    ind = gpDataIndices(model, i)
	    ll =  ll - (sum(log(diag(model$Lm[[i]])))
		      + 0.5*(-(length(ind) - model$k)*log(model$beta)
		      +(length(ind)*log(2*pi)+sum(log(model$diagD[[i]])))))
	    ll = ll - 0.5*model$beta*(t(model$scaledM[[i]])%*%model$scaledM[[i]]
		    - t(model$bet[[i]])%*%model$bet[[i]])
	  }
	}
      } else
	stop("Variable length Beta not implemented yet.")
    }
  } else if (model$approx == "pitc") {
    ## Partially independent training conditional.
    if ((!"isSpherical" %in% names(model)) || model$isSpherical) {
      if (length(model$beta)==1) {
	ll = model$d*(model$logDetA-model$logDetK_uu + model$k*log(model$beta))
	## Loop through the blocks computing each part to be added.
	K_ufDinvm = matrix(0, model$k, model$d)
	Dinvm = list()
	for (i in 1:length(model$blockEnd)) {
	  ind = gpBlockIndices(model, i)
	  Dinvm[[i]] = model$Dinv[[i]]%*%model$m[ind, ,drop=FALSE]
	  K_ufDinvm = K_ufDinvm + model$K_uf[, ind,drop=FALSE]%*%Dinvm[[i]]
	}
	ll = ll - model$beta*sum((model$Ainv%*%K_ufDinvm) * K_ufDinvm)
	
	for (i in 1:length(model$blockEnd)) {
	  ind = gpBlockIndices(model, i)
	  ll = ll + model$d*(model$logDetD[i] - length(ind)*log(model$beta))
		  + model$beta*sum(Dinvm[[i]] * model$m[ind, ,drop=FALSE])
	}
	ll = -0.5*ll
	ll = ll - 0.5*model$N*model$d*log(2*pi)
      } else
	stop("Variable Length Beta not implemented yet.")
    } else {
      if (length(model$beta)==1) {
	ll = 0
	Dinvm = matrix(0, model$blockEnd, model$d)
	Dinvm = lapply(split(Dinvm,row(Dinvm)), split, 1:model$d)
	for (j in 1:model$d) {
	  ll = ll + model$logDetA[j]-model$logDetK_uu + model$k*log(model$beta)
	  ## Loop through the blocks computing each part to be added.
	  K_ufDinvm = matrix(0, model$k, 1)
	  for (i in 1:length(model$blockEnd)) {
	    ind = gpDataIndices(model, j, i)
	    Dinvm[[i]][[j]] = model$Dinv[[i]][[j]]%*%model$m[ind, j,drop=FALSE]
	    K_ufDinvm = K_ufDinvm + model$K_uf[, ind]%*%Dinvm[[i]][[j]]
	  }
	  ll = ll - model$beta*sum((model$Ainv[[i]]%*%K_ufDinvm) * K_ufDinvm)

	  for (i in 1:length(model$blockEnd)) {
	    ind = gpDataIndices(model, j, i)
	    ll = ll + model$logDetD[i, j] - length(ind)*log(model$beta)
		    + model$beta*sum(Dinvm[[i]][[j]] * model$m[ind, j,drop=FALSE])
	    ll = ll + length(ind)*log(2*pi)
	  }
	}
	ll = -0.5*ll
      } else
	stop("Variable Length Beta not implemented yet.")
    }
  }

  if (model$learnScales)
    ll = ll - sum(log(model$scale))

  ll = ll - model$d * model$N/2 * log(2*pi)
    
  return (ll)
}