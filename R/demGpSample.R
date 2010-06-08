## DEMGPSAMPLE Simple demonstration of sampling from a covariance function.

bw = FALSE
x = as.matrix(seq(-1,1,length=25)) #linspace(-1, 1, 25)'
kern = kernCreate(x, 'rbf')
kern$inverseWidth = 10

par(pty="s") #figure(1)
par(mar=c(4, 4, 2, 6)) ## increase right margin for legend
K = kernCompute(kern, x)
colormap=heat.colors(256)
if (bw)
  colormap=gray.colors(256)

image(t(K), xlab='m', ylab='n')
#colorbar
image.plot(legend.only = TRUE, nlevel = 256, zlim = c(min(K),max(K)), col = colormap)

# t = [];
# t = [t xlabel('n')];
# t = [t ylabel('m')];
# set(gca, 'fontname', 'times')
# set(gca, 'fontsize', 24)
if (exists('printDiagram') && exists('var') && exists(printDiagram))
  printPlot('gpCovariance', '../tex/diagrams/', '../html') ## !!!

## Need to take the real part of the sample, as the kernel is numerically less than full rank.
f = t(Re(gaussSamp(K, 1)))

dev.new() #figure(2) 
a = plot(f, xlab='n', ylab= bquote(f[n]))

# t = [t text(12.5, -.5, 'n')];
# t = [t text(-5, .5, 'f_n')];
# set(t, 'fontname', 'times')
# set(t, 'fontsize', 24)
# set(t, 'fontangle', 'italic')
# set(gca, 'fontname', 'times')
# set(gca, 'fontsize', 18)
# set(a,'markersize', 20)
# set(a, 'linewidth', 2)
# zeroAxes(gca)
if (exists('printDiagram') && exists('var') && exists(printDiagram))
  printPlot('gpSample', '../tex/diagrams/', '../html') ## !!!

# save(K, x, f, file='demGpSample.RData')