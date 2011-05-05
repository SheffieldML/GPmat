hmcDefaultOptions <- function() {
  return (list(tau=10, epsilon=0.05, maxit=1000, threshold=-36))
}


HMCsample <- function (x, fn, grad, options, ...) {
  L <- options$maxit

  dim <- dim(x)
  xhist <- matrix(0, nrow=L, ncol=dim[2])
  Ehist <- rep(0, L)

  Tau <- options$tau
  epsilon <- options$epsilon
  threshold <- options$threshold
  
  g <- grad(x, ...)
  E <- fn(x, ...)

  for (l in seq(L)) {
    p <- matrix(rnorm(prod(dim)), nrow=dim[1], ncol=dim[2])
    H <- p %*% t(p) / 2 + E    # evaluate H(x, p)
    xnew <- x ; gnew <- g
    for (tau in seq(Tau)) {
      p <- p - epsilon * gnew / 2 ; # make half-step in p 
      xnew <- xnew + epsilon * p ;  # make step in x 
      gnew <- grad(xnew, ...);      # find new gradient 
      p <- p - epsilon * gnew / 2 ; # make half-step in p 
    }
    xnew[xnew < threshold] <- threshold
    Enew <- fn(xnew, ...) ;         # find new value of H
    Hnew <- p %*% t(p) / 2 + Enew ; 
    dH <- Hnew - H ;                # Decide whether to accept 
    cat(sprintf('Step %d, threshold: %.4f\n', l, exp(-dH)));
    if ( runif(1) < exp(-dH) ) {
      cat('accepted.\n');
      cat(x);
      cat('\n');
      g <- gnew ; x <- xnew ; E <- Enew ;
    }
    xhist[l,] <- x;
    Ehist[l] <- E;
  }
  return (list(x=x, E=E, xhist=xhist, Ehist=Ehist))
}


modelSample <- function (model, options, ...) {
  if (is.GPModel(model)) {
    haveGPModel <- TRUE
    model <- modelStruct(model)
  }
  else
    haveGPModel <- FALSE

  params <- modelExtractParam(model)
  samples <- HMCsample(params, modelObjective, modelGradient, options, model)

  model <- modelExpandParam(model, samples$x)
  model$llscore <- samples$E

  if (haveGPModel)
    return (list(model=new("GPModel", model), samples=samples))
  else
    return (list(model=model, samples=samples))
}
