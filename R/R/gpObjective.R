gpObjective <- function(params, model) {

  model = gpExpandParam(model, params)
  f = - gpLogLikelihood(model)

  return (f)
}
