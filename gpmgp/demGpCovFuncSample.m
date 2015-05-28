function demGpCovFuncSample

% DEMGPCOVFUNCSAMPLE Sample from some different covariance functions.
% FORMAT
% DESC samples from some different covariance functions with different
% kernel parameters.
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008
%
% SEEALSO : kernCreate
  
% GP

  global printDiagram
  randn('seed', 1e5)
  rand('seed', 1e5)
  
  x = linspace(-2, 2, 200)';
  kern = kernCreate(x, 'rbf');
  
  numSamps = 10;
  
  figure(1);
  kern.inverseWidth = 1/(0.3^2);
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);
  
  figure(2)
  kern.inverseWidth = 1;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);
  
  figure(3)
  kern.inverseWidth = 1/(0.3^2);
  kern.variance = 4;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);
  
  figure(4)
  kern = kernCreate(x, 'bias');
  kern.variance = 4;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);

  
  figure(5)
  kern = kernCreate(x, 'lin');
  kern.variance = 16;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);

  figure(6)
  kern = kernCreate(x, 'poly');
  kern.variance = 1;
  kern.biasVariance = 1;
  kern.weightVariance = 1;
  kern.degree = 5;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);

  figure(7)
  kern = kernCreate(x, 'mlp');
  kern.variance = sqrt(pi)/2;
  kern.biasVariance = 100;
  kern.weightVariance = 100;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);
  
  
  figure(8)
  kern = kernCreate(x, 'mlp');
  kern.variance = sqrt(pi/2);
  kern.biasVariance = 0;
  kern.weightVariance = 100;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);
  
  
  figure(9)
  kern = kernCreate(x, {'rbf', 'bias', 'white'});
  kern.comp{1}.variance = 1;
  kern.comp{1}.inverseWidth = 10;
  kern.comp{2}.variance = 1;
  kern.comp{3}.variance = 0.01;
  K = kernCompute(kern, x);
  f = real(gsamp(zeros(1, size(x, 1)), K, numSamps))';
  a = plot(x, f);
  prepPlot(a);

end


function prepPlot(a)

  global printDiagram
  set(gca, 'xlim', [-2 2])
  set(gca, 'ylim', [-4 4])
  set(gca, 'fontname', 'times')
  set(gca, 'fontsize', 18)
  set(a,'markersize', 10)
  set(a, 'linewidth', 2)
  zeroAxes(gca, 0.025, 18, 'times');
  fig = gcf;
  if printDiagram
    printPlot(['demGpCovFuncSample' num2str(fig)], '../tex/diagrams', '../html');
  end
end
