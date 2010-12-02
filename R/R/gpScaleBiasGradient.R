gpScaleBiasGradient <- function(model) {
  g = list()
  if (model$learnScales) {
    ## 'drop' converts row matrix to column vector by default.
    g = 1/model$scale * drop(model$innerProducts-1)
    fhandle <- get(model$scaleTransform$func, mode="function")
    g = g * fhandle(model$scale, "gradfact")
  }

  return (g)
}