% DEMGPSAMPLE Simple demonstration of sampling from a covariance function.

% GP

randn('seed', 1e5)
rand('seed', 1e5)

bw = false;
x = linspace(-1, 1, 25)';
kern = kernCreate(x, 'rbf');
kern.inverseWidth = 10;

figure(1)
clf
K = kernCompute(kern, x);
imagesc(K);
if bw
  colormap gray
else
  colormap jet
end  
colorbar
t = [];
t = [t xlabel('n')];
t = [t ylabel('m')];
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 24)
if exist('printDiagram', 'var') && printDiagram
  printPlot('gpCovariance', '../tex/diagrams/', '../html');
end
% need to take the real part of the sample as the kernel is numerically less than full rank 
f = real(gsamp(zeros(1, size(x, 1)), K, 1))';

figure(2) 
clf
a = plot(f, 'k.');
t = [t text(12.5, -.5, 'n')];
t = [t text(-5, .5, 'f_n')];
set(t, 'fontname', 'times')
set(t, 'fontsize', 24)
set(t, 'fontangle', 'italic')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
set(a,'markersize', 20)
set(a, 'linewidth', 2)
zeroAxes(gca)
if exist('printDiagram', 'var') && printDiagram
  printPlot('gpSample', '../tex/diagrams/', '../html');
end
% if bw
%   set(a, 'color', [0 0 0])
%   print('-deps', '../tex/diagrams/gpSamplebw.eps');
% else
%   print('-depsc', '../tex/diagrams/gpSample.eps');
% end

save demGpSample K x f
