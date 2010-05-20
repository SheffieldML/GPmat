gpLogLikeGradients <- function(model, X, M, X_u, gX_u.return=FALSE, gX.return=FALSE,
  g_beta.return=FALSE) {

  if (nargs() < 4) {
    if ("X_u" %in% names(model))
      X_u = model$X_u
    else
      X_u = list()

    if (nargs()< 3 && (!"S" %in% names(model)))
      M = model$m

    if (nargs() < 2)
      X = model$X
  }

  gX_u = list()
  gX = list()

  g_scaleBias = gpScaleBiasGradient(model)
  g_meanFunc = list()
  if ("meanFunction" %in% names(model) && length(model$meanFunction)>0)
    g_meanFunc = gpMeanFunctionGradient(model)    

  if (model$approx == "ftc") {
    ## Full training conditional.
    if (gX_u.return && gX.return) {
#       ## Prepare to Compute Gradients with respect to X
#       gKX = kernGradX(model$kern, X, X) ## !!! 
#       gKX = gKX*2
#       dgKX = kernDiagGradX(model$kern, X) ## !!!
#       for (i in 1:model$N)
#         gKX[i, , i] = dgKX[i, ]
# 
#       gX = matrix(0, model$N, model$q)
    }

    ## Gradients of Kernel Parameters
    g_param = matrix(0, 1, model$kern$nParams)
    g_beta = list()
    if ("beta" %in% names(model))
      g_beta = 0

    ## For very high D, we use the matrix S which is M*M'
    if ("S" %in% names(model)) {
#       gK = localSCovarianceGradients(model) ## !!!
#       if (gX_u.return && gX.return) {
#         ## Compute Gradients with respect to X
#         counter = 0
#         for (i in 1:model$N) {
#           counter = counter + 1
#           for (j in 1:model$q)
#             gX[i, j] = gX[i, j] + gKX[, j, i] %*% gK[, counter]
#         }
#       }
#       ## Compute Gradients of Kernel Parameters
#       g_param = g_param + kernGradient(model$kern, X, gK) ## !!!
    } else {
      for (k in 1:model$d) {
        gK = localCovarianceGradients(model, M[, k], k)
        if (gX_u.return && gX.return) {
#           ## Compute Gradients with respect to X
#           ind = gpDataIndices(model, k)
#           counter = 0
#           for (i in ind) {
#             counter = counter + 1
#             for (j in 1:model$q)
#               gX[i, j] = gX[i, j] + gKX[ind, j, i]%*%gK[, counter]
#           }
        }
        ## Compute Gradients of Kernel Parameters
        if (model$isMissingData)
#           g_param = g_param
# 		+ kernGradient(model$kern, X[model$indexPresent[[k]], ], gK) ## !!!
        else
          g_param = g_param + kernGradient(model$kern, X, gK)
      }

      if ("beta" %in% names(model) && model$optimiseBeta) {
#         if (dim(model$beta)[1] == 1)
#           g_beta = g_beta + sum(diag(gK))
#         else if (dim(model$beta)[2]==1 && dim(model$beta)[1]==model$N)
#           g_beta = g_beta + diag(gK)
#         else if (dim(model$beta)[2]==model$d && dim(model$beta)[1]==model$N)
#           g_beta[, k] = diag(gK)
#         else
#           stop("Unusual dimensions for model.beta.")
      }
    }
  } else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
#     ## Sparse approximations.
#     ## [gK_u, gK_uf, gK_star, g_beta] = gpCovGrads(model, M) ## !!!
#     gK = gpCovGrads(model, M) ## !!!
#     gK_u=gk$gk_u; gK_uf=gK$gK_uf; gK_star=gk$gK_star; g_beta=gK$g_beta
# 
#     ## Compute Gradients of Kernel Parameters
#     gParam_u = kernGradient(model$kern, X_u, gK_u)
#     gParam_uf = kernGradient(model$kern, X_u, X, gK_uf)
# 
#     g_param = gParam_u + gParam_uf
# 
#     ## Compute Gradients with respect to X_u
#     gKX = kernGradX(model$kern, X_u, X_u) ## !!!
# 
#     ## The 2 accounts for the fact that covGrad is symmetric
#     gKX = gKX*2
#     dgKX = kernDiagGradX(model$kern, X_u) ## !!!
#     for (i in 1:model$k)
#       gKX[i, , i] = dgKX[i, ]
#     
#     if (!model$fixInducing || gX_u.return || gX.return || g_beta.return) { #nargout > 1
#       ## Allocate space for gX_u
#       gX_u = matrix(0, model$k, model$q)
#       ## Compute portion associated with gK_u
#       for (i in 1:model$k) {
#         for (j in 1:model$q)
#           gX_u[i, j] = gKX[, j, i]%*%gK_u[, i]
#       }
#       
#       ## Compute portion associated with gK_uf
#       gKX_uf = kernGradX(model$kern, X_u, X) ## !!!
#       for (i in 1:model$k) {
#         for (j in 1:model$q)
#           gX_u[i, j] = gX_u[i, j] + gKX_uf[, j, i]%*%gK_uf[i, ]
#       }
#     }
# 
#     if (gX_u.return && gX.return) { #nargout > 2
#       ## Compute gradients with respect to X
#       ## Allocate space for gX
#       gX = matrix(0, model$N, model$q)
#       
#       ## this needs to be recomputed so that it is wrt X not X_u
#       gKX_uf = kernGradX(model$kern, X, X_u) ## !!!
#       
#       for (i in 1:model$N) {
#         for (j in 1:model$q)
#           gX[i, j] = gKX_uf[, j, i]%*%gK_uf[, i]
#       }
#     }
  } else
    stop("Unknown model approximation.")

  if (model$approx == "ftc") {
    ## Full training conditional. Nothing required here.
  } else if (model$approx == "dtc") {
    ## Deterministic training conditional.  
  } else if (model$approx %in% c("fitc","dtcvar")) {
#     ## Fully independent training conditional.
#     ## Variational sparse approximation.
#     
#     if (gX_u.return && gX.return) { #nargout > 2
#       ## deal with diagonal term's effect on X gradients.
#       gKXdiag = kernDiagGradX(model$kern, X) ## !!!
#       for (i in 1:model$N)
#         gX[i, ] = gX[i, ] + gKXdiag[i, ]%*%gK_star[i]
#     }
#     
#     ## deal with diagonal term's affect on kernel parameters.
#     g_param = g_param + kernDiagGradient(model$kern, X, gK_star) ## !!!
  } else if (model$approx == "pitc") {
#     ## Partially independent training conditional.    
#     if (gX_u.return && gX.return) { #nargout > 2
#       ## deal with block diagonal term's effect on X gradients.
#       startVal = 1
#       for (i in 1:length(model$blockEnd)) {
#         endVal = model$blockEnd[i]
#         ind = startVal:endVal
#         gKXblock = kernGradX(model$kern, X[ind, ], X[ind, ]) ## !!!
#         ## The 2 accounts for the fact that covGrad is symmetric
#         gKXblock = gKXblock*2
# 
#         ## fix diagonal
#         dgKXblock = kernDiagGradX(model$kern, X[ind, ]) ## !!!
#         for (j in 1:length(ind))
#           gKXblock[j, , j] = dgKXblock[j, ]
#         
#         for (j in ind) {
#           for (k in 1:model$q) {
#             subInd = j - startVal + 1
#             gX[j, k] = gX[j, k] + gKXblock[, k, subInd]%*%gK_star[[i]][, subInd]
#           }
#         }
#         startVal = endVal + 1
#       }
#     }
#     ## deal with block diagonal's effect on kernel parameters.
#     for (i in 1:length(model$blockEnd)) {
#       ind = gpBlockIndices(model, i)
#       g_param = g_param + kernGradient(model$kern, X[ind, ], gK_star[[i]])
#     }
  } else
    stop("Unrecognised model approximation")
  
  if (!(gX_u.return && gX.return && g_beta.return)) { #if (nargout < 4)
    if ((!"optimiseBeta" %in% names(model) && model$approx!="ftc") || model$optimiseBeta)
      ## append beta gradient to end of parameters
      gParam = c(g_param, g_meanFunc, g_scaleBias, g_beta)
    else
      gParam = c(g_param, g_meanFunc, g_scaleBias)
  } else
    gParam = c(g_param, g_meanFunc g_scaleBias)

  ## if there is only one output argument, pack gX_u and gParam into it.
  if (!(gX_u.return || gX.return || g_beta.return)) { #(nargout == 1)
    gParam = c(gX_u, gParam)

  return (gParam)
} # function




localCovarianceGradients <- function(model, y, dimension) {
  if (!"isSpherical" %in% names(model) || model$isSpherical) {
    invKy = model$invK_uu %*% y
    gK = -model$invK_uu + invKy%*%invKy
  } else {
    if (model$isMissingData)
      m = y[model$indexPresent[[dimension]]]
    else
      m = y

    invKy = model$invK_uu[[dimension]] %*% m
    gK = -model$invK_uu[[dimension]] + invKy%*%invKy
  }
  gK = gK * .5
  return (gK)
}

function gK = localSCovarianceGradients(model) {
  gK = -model$d%*%model$invK_uu + model$invK_uu%*%model$S%*%model$invK_uu
  gK = gK * .5
  return (gK)
}