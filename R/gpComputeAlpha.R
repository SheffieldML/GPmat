gpComputeAlpha <- function(model, m) {

  if (nargs() < 2)
    m = model$m

  if (model$approx == "ftc") {
    model$alpha = matrix(0,model$N,model$d) #zeros(model.N, model.d)
    if (!"isSpherical" %in% names(model) || model$isSpherical)
      model$alpha = model$invK_uu %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[ind, i] = model$invK_uu[[i]] %*% m[ind, i]
      }
    }
  }
  else if (model$approx %in% c("dtc","dtcvar")) {
    model$alpha = matrix(0, model$k, model$d)
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]] %*% model$K_uf[,ind] %*% m[ind,i]
      }
    }
  }
  else if (model$approx == "fitc") {
    model$alpha = matrix(0, model$k, model$d)
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% model$Dinv %*% m
    else {
      for (i in 1:model$d) {
	ind = gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]]
	    %*% model$K_uf[,ind] %*% model$Dinv[[i]] %*% m[ind,i]
      }
    }
  else if (model$approx == "pitc") {
    model$alpha = matrix(0, model$k, model$d)
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      for (i in seq(along=model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	model$alpha = model$alpha + model$Ainv%*%model$K_uf[,ind]%*%model$Dinv[[i]]%*%m[ind,]
      }
    else {
      for (i in seq(along=model.blockEnd)) {
	for (j in 1:model$d) {
	  ind = gpDataIndices(model, j, i)
	  model$alpha[,j] = model$alpha[,j] + model$Ainv[[j]]%*%model$K_uf[,ind]%*%
	      model$Dinv[[i,j]]%*%m[ind,j]
	}
      }
    }
  }

  return(model)
}