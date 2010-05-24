gpOptimise <- function(model, display=TRUE, iters=2000, gradcheck=FALSE) {

#   if (nargs() < 4) {
#     gradcheck = FALSE
#     if (nargs() < 3) {
#       iters = 2000
#       if (nargs() < 2)
# 	display = TRUE
#     }
#   }

  params = gpExtractParam(model)
  ## options = list(maxit=3000, ln=c(0,2), xtol=1e-4, fnTol=1e-4, optimiser="SCG",
  ##	gradcheck=FALSE, display=TRUE)
  options = optimiDefaultOptions() #optOptions
  if (display) {
    options$display = TRUE
    if ((length(params) <= 100) && gradcheck)
      options$gradcheck = TRUE
  }
  options$maxit = iters

  if ("optimiser" %in% names(model))
    optim = get(paste(model$optimiser, "optim", sep=""), mode="function")
  else
    optim = get("CGoptim", mode="function")

  fn = get('gpObjective', mode="function")
  grad = get('gpGradient', mode="function")

#   strcmp(func2str(optim), 'optimiMinimize')
#     ## Carl Rasmussen's minimize function 
#     params = optim('gpObjectiveGradient', params, options, model);
#   else
  ## R version of NETLAB function
  ## only necessary parts from gpGradient implemented!
  params = optim(params, fn, grad, options, model) ## !!!

  model = gpExpandParam(model, params$xmin)

  return (model)
}