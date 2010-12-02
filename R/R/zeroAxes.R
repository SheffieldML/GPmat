zeroAxes <- function(xy=NULL, col='blue') {
  abline(v = 0, h = 0, col=col, lwd=.5) ## Axes
#   scx = (max(xy[,1])-min(xy[,1]))/100
#   scy = (max(xy[,2])-min(xy[,2]))/100
#   lines(scx*cos(seq(0, 2*pi, l=100)), (scx)*sin(seq(0, 2*pi, l=100)),col=col,lwd=.5)
}