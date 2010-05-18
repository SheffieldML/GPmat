gpScaleBiasGradient <- function(model) {
  g = list()
  if (model$learnScales) {
    g = 1/model$scale * (model$innerProducts-1)
    fhandle <- get(paste(model$scaleTransform, "Transform", sep=""), mode="function")
    g = g * fhandle(model$scale, "gradfact")
  }

  return (g)
}