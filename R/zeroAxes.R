zeroAxes <- function(xy, col='blue') {
  abline(v = 0, h = 0, col=col, lwd=.5) ## axes and co-centric circles
  int = (max(xy[,1])-min(xy[,1]))/100
  lines(int*cos(seq(0, 2*pi, l=100)), int*sin(seq(0, 2*pi, l=100)),col=col,lwd=.5)
}