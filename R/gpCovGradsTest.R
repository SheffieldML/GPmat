gpCovGradsTest <- function(model) {
   
  changeVal = 1e-6
  if (model$approx == "ftc") {
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    diffsK_uu = diag(0, dim(model$K_uu)[1])
    diffsK_uf = matrix(0, dim(model$K_uf)[1], dim(model$K_uf)[2])
    for (i in 1:dim(model$K_uu)[1]) {
      for (j in 1:i) {
	origK = model$K_uu[i, j]
	model$K_uu[i, j] = origK + changeVal
	model$K_uu[j, i] = model$K_uu[i, j]
	invK = .jitCholInv(model$K_uu, silent=TRUE)
	model$invK_uu = invK$invM
	model$logDetK_uu = 2*sum(log(diag(invK$chol)))
	model = gpUpdateAD(model)
	objPlus = gpLogLikelihood(model)
	model$K_uu[i, j] = origK - changeVal
	model$K_uu[j, i] = model$K_uu[i, j]
	## re-compute gradients at same location + offset
	invK = .jitCholInv(model$K_uu, silent=TRUE)
	model$invK_uu = invK$invM
	model$logDetK_uu = 2*sum(log(diag(invK$chol)))
	model = gpUpdateAD(model)
	objMinus = gpLogLikelihood(model)
	## compute finite differences
	diffsK_uu[i, j] = (objPlus - objMinus)/(2*changeVal)
	diffsK_uu[j, i] = diffsK_uu[i, j]
	model$K_uu[i, j] = origK
	model$K_uu[j, i] = origK
	invK = .jitCholInv(model$K_uu, silent=TRUE)
	model$invK_uu = invK$invM
	model$logDetK_uu = 2*sum(log(diag(invK$chol)))
	model = gpUpdateAD(model)
      }
    }

    for (i in 1:dim(model$K_uf)[1]) {
      for (j in 1:dim(model$K_uf)[2]) {
	origK = model$K_uf[i, j]
	model$K_uf[i, j] = origK + changeVal
	model = gpUpdateAD(model)
	objPlus = gpLogLikelihood(model)
	model$K_uf[i, j] = origK - changeVal
	model = gpUpdateAD(model)
	objMinus = gpLogLikelihood(model)
	diffsK_uf[i, j] = (objPlus - objMinus)/(2*changeVal)
	model$K_uf[i, j] = origK
	model = gpUpdateAD(model)
      }
    }

    temp = gpCovGrads(model, model$m)
    gK_uu = temp$gK_uu; gK_uf = temp$gK_uf; g_Lambda = temp$g_Lambda
    gK_uuMaxDiff = max(abs(2*(gK_uu-diag(diag(gK_uu))) + diag(diag(gK_uu)) - diffsK_uu))
    gK_ufMaxDiff = max(abs(gK_uf - diffsK_uf))

    print(paste("K_uu grad max diff ", round(gK_uuMaxDiff,4), sep=""))
    if (gK_uuMaxDiff > 1e-4)
      print(2*(gK_uu-diag(diag(gK_uu))) + diag(diag(gK_uu)) - diffsK_uu)

    print(paste("K_uf grad max diff ", round(gK_ufMaxDiff,4), sep=""))
    if (gK_ufMaxDiff > 1e-4)
      print(gK_uf - diffsK_uf)
  }

  return (model)
}
