gpGradient <- function(params, model) {

  model = gpExpandParam(model, params)
  ## only necessary parts from gpLogLikeGradients implemented!
  g = - gpLogLikeGradients(model)

  return (g)
}
