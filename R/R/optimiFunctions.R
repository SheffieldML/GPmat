optimiDefaultConstraint <- function (constraint) {
  if ( constraint == "positive" ) {
    return (list(func="expTransform", hasArgs=FALSE))
  } else if ( constraint == "zeroone" ) {
    return (list(func="sigmoidTransform", hasArgs=FALSE))
  } else if ( constraint == "bounded" ) {
    return (list(func="boundedTransform", hasArgs=TRUE))
  }
}



## optimiFdefaultOptions <- function () {
##   options <- list(maxit=100, tolf=1e-4, tol=1e-4, display=TRUE)
##   return (options)
## }



optimiDefaultOptions <- function() {
  ## options: trace, maximum iteration, fnscale, reltol, default optimiser
  ## return (list(trace=TRUE, maxit=1000, fnscale=1e1, reltol=1e-4, optimiser="CG", gradcheck=FALSE, hessian=FALSE))
  return (list(maxit=3000, ln=c(0,2), xtol=1e-4, fnTol=1e-4, optimiser="SCG", gradcheck=FALSE, display=TRUE))
}



modelOptimise <- function (model, options, ...) {
  if (is.GPModel(model)) {
    haveGPModel <- TRUE
    model <- modelStruct(model)
  }
  else
    haveGPModel <- FALSE

  funcName <- paste(model$type, "Optimise", sep="")
  if ( exists(funcName, mode="function") ) {
    func <- get(funcName, mode="function")
    model <- func(model, options, ...)
  } else {
    if ( "optimiser" %in% names(options) ) {
      funcName <- paste(options$optimiser, "optim", sep="")
    } else {
      funcName <- "CGoptim"
    }
    optFunc <- get(funcName, mode="function")

    params <- modelExtractParam(model)
    newParams <- optFunc(params, modelObjective, modelGradient, options, model)
    
    model <- modelExpandParam(model, newParams$xmin)
    model$llscore <- newParams$objective
  }

  if (haveGPModel)
    return (new("GPModel", model))
  else
    return (model)
}



.fn_line <- function (linemin, fun, para0, direction, ...) {
  ## y = fn (x)
  func <- function(x, ...) fun(x, ...)
  
  para <- para0 + linemin * direction
  
  ans <- func(para, ...)
  
  return (ans)
}


CGoptim <- function (x, fn, grad, options, ...) {
  ## option[1] : number of iterations
  ## option[2] : interval for the line search
  ## option[3] : tolerence for x to terminate the loop
  ## option[4] : tolerence for fn to terminate the loop
  ## option$display : option of showing the details of optimisaton

  ## y = fn (x)
  func <- function(x, ...) fn(x, ...)
      
  ## gradient function = gr (x)
  gradfunc <- function(x, ...) grad(x, ...)
	
  fn_new <- func(x, ...)
  #if ( display ) 
  #  cat ("fn0 :",fn_new, "\n")

  grad_new <- gradfunc(x, ...)
  #if ( options$display )  
  #  cat ("grad0 :",grad_new, "\n\n")
	
  direction <- -grad_new
  lnSchFail <- FALSE
  for ( ind in 1:options[[1]] ) {
	
    x_old <- x
    fn_old <- fn_new
    grad_old <- grad_new
    
    grad2 <- crossprod(grad_old)	
    
    if ( grad2 == 0 ) {
      objective <- fn_new
      xmin <- x
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }

    dnorm <- sqrt(sum(direction*direction))
    line_dir <- direction / dnorm
    ## cat ("\n line_dir :", line_dir, "\n\n")
    ow <- options("warn")
    options(warn=2)
    lnSch <- try( optimize(.fn_line, options[[2]], para0=x_old, direction=line_dir, fun=fn, ...), silent=TRUE )
    options(ow)

    if ( is.list(lnSch) ) {
      x <- x_old + lnSch$minimum * line_dir		
      fn_new <- lnSch$objective
      fnmin <- min(fn_old, fn_new)
      if ( fnmin==fn_new ) {
        xnmin <- x
      } else {
        xnmin <- x_old
      }
    } else {
      warning("Line search failed! \n")
      x <- xnmin
      fn_new <- fnmin
      lnSchFail <- TRUE

      xmin <- xnmin
      objective <- fnmin
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }
    
    if ( max(abs(x-x_old))<options[[3]] & max(abs(fn_new-fn_old))<options[[4]] ) {
      xmin <- x
      objective <- fn_new
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }
    
    grad_new <- gradfunc(x, ...)
    
    eta <- ( t(grad_new-grad_old) %*% grad_new ) / grad2
    direction <- direction * eta - grad_new

    if ( options$display )
      cat(ind, "-th objective = :", fn_new, "\t max xi: ", max(abs(x-x_old)), "\n")
  }

  warning("Maximum iteration reached! \n")
  xmin <- x
  objective <- fn_new
  
  ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
  ##if ( options$display ) 
  ##  cat("\n Optimisation ends! \n")
  return (ans)
}



SCGoptim <- function (x, fn, grad, options, ...) {
  ## options = list(maxit, ln, xtol, fnTol, optimiser="SCG", gradcheck=FALSE)
  ##cat ("\n SCG Optimisation begins! \n")

  if ( "maxit" %in% names(options) && !is.na(options$maxit) ) {
    niters <- options$maxit
  } else {
    niters <- 100
  }

  display <- options$display
  gradcheck <- options$gradcheck
  
  ## y = fn (x)
  func <- function(x, ...) fn(x, ...)
      
  ## gradient function = gr (x)
  gradfunc <- function(x, ...) grad(x, ...)

  nparams <- length(x)

  sigma0 <- 1e-4
  fold <- func(x, ...)
  fnow <- fold
  gradnew <- gradfunc(x, ...)
  gradold <- gradnew
  d <- -gradnew
  success <- 1
  nsuccess <- 0
  beta <- 1
  betamin <- 1e-15
  betamax <- 1e100
  eps <- 2.2204e-16
  j <- 1

  while ( j<=niters ) {
    if ( success == 1 ) {
      mu <- crossprod(d, gradnew)
      if ( mu>=0 ) {
        d <- -gradnew
        mu <- crossprod(d, gradnew)
      }
      kappa <- crossprod(d, d)
      if ( kappa<eps ) {
        xmin <- x
        objective <- fnow
        ans <- list(xmin=xmin, objective=objective)
        return (ans)
      }
      sigma <- (sigma0/sqrt(kappa))[1]
      xplus <- x+sigma*d
      gplus <- gradfunc(xplus, ...)
      theta <- crossprod(d, gplus-gradnew)/sigma
    }

    delta <- theta + beta*kappa
    if ( delta<=0 ) {
      delta <- beta*kappa
      beta <- beta - theta/kappa
    }
    alpha <- (-mu/delta)[1]

    xnew <- x+alpha*d
    ow <- options("warn")
    options(warn=2)
    fnew <- try( func(xnew, ...), silent=TRUE )
    options(ow)
    if ( !is.finite(fnew) ) fi <- 1
    while ( !is.finite(fnew) ) {
      alpha <- alpha/2
      xnew <- x+alpha*d
      ow <- options("warn")
      options(warn=2)
      fnew <- try( func(xnew, ...), silent=TRUE )
      options(ow)
      fi <- fi+1
      if ( is.finite(fnew) ) {
        fi <- 0
      }
    }

    Delta <- 2*(fnew-fold)/(alpha*mu)
    if ( Delta>=0 ) {
      success <- 1
      nsuccess <- nsuccess+1
      x <- xnew
      fnow <- fnew
    } else {
      success <- 0
      fnow <- fold
    }
    if (display)
      cat("Cycle ", j, "Error ", round(fnow, digits=4), "Scale ", beta, "\n")

    if ( success == 1 )
      if ( (max(abs(alpha*d))<options$xtol) & (max(abs(fnew-fold))<options$fnTol) ) {
        xmin <- x
        objective <- fnew
        ans <- list(xmin=xmin, objective=objective)
        return (ans)
      } else {
        fold <- fnew
        gradold <- gradnew
        gradnew <- gradfunc(x, ...)
        if ( crossprod(gradnew, gradnew) == 0 ) {
          xmin <- x
          objective <- fnew
          ans <- list(xmin=xmin, objective=objective)
          return (ans)
        }
      }

    if ( Delta < 0.25 ) 
      beta <- min(2*beta, betamax)

    if ( Delta > 0.75 )
      beta <- max(0.5*beta, betamin)

    if ( nsuccess == nparams ) {
      d <- -gradnew
      nsuccess <- 0
    } else {
      if ( success==1 ) {
        gamma <- (crossprod(gradold-gradnew, gradnew)/mu)[1]
        d <- gamma*d-gradnew
      }
    }

    j <- j+1
  }

  xmin <- x
  objective <- fold
  ans <- list(xmin=xmin, objective=objective)
  warning("Maximum number of iterations has been exceeded.\n")
  return (ans)
}
