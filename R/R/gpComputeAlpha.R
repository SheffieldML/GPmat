gpComputeAlpha <- function(model, m) {

  if (nargs() < 2)
    m = model$m

  model$alpha = matrix(0, model$k, model$d)
  if (model$approx == "ftc") {
# browser()	## m = y-mean(y)
    if (!"isSpherical" %in% names(model) || model$isSpherical)
      model$alpha = model$invK_uu %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[ind, i] = model$invK_uu[[i]] %*% m[ind, i,drop=FALSE]
      }
    }
  }
  else if (model$approx %in% c("dtc","dtcvar")) {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]] %*% model$K_uf[,ind,drop=FALSE] %*% m[ind,i,drop=FALSE]
      }
    }
  }
  else if (model$approx == "fitc") {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% model$Dinv %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]] %*%
			    model$K_uf[,ind,drop=FALSE] %*% model$Dinv[[i]] %*% m[ind,i,drop=FALSE]
      }
    }
  } else if (model$approx == "pitc") {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      for (i in seq(along=model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	model$alpha = model$alpha + model$Ainv%*%model$K_uf[,ind,drop=FALSE]%*%model$Dinv[[i]]%*%m[ind, ,drop=FALSE]
      }
    else {
      for (i in seq(along=model$blockEnd)) {
	for (j in 1:model$d) {
	  ind = gpDataIndices(model, j, i)
	  model$alpha[,j] = model$alpha[,j,drop=FALSE] + model$Ainv[[j]]%*%model$K_uf[,ind,drop=FALSE]%*%
	      model$Dinv[[i]][[j]]%*%m[ind,j,drop=FALSE]
	}
      }
    }
  }

  return(model)
}
