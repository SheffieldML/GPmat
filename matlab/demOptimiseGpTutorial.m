% DEMOPTIMISEGPTUTORIAL Shows that there is an optimum for the covariance function length scale.
% DESC shows that by varying the length scale an artificial data
% set has different likelihoods, yet there is an optimum for which
% the likelihood is maximised.

% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% GP

randn('seed', 1e5);
rand('seed', 1e5);

fillColor = [0.7 0.7 0.7];
markerSize = 40;
markerWidth = 4;
markerType = 'k.';
lineWidth = 4;


x = linspace(-1, 1, 6)';
trueKern = kernCreate(x, {'rbf', 'white'});
kern.comp{2}.variance = 0.001;
K = kernCompute(trueKern, x);
y = gsamp(zeros(1, 6), K, 1)';

xtest = linspace(-2, 2, 200)';
kern = trueKern;

lengthScale = [0.05 0.2 0.5 1 2 10];
counter = 0;

figure(1)
p = plot(x, y, markerType);
set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 24)
set(gca, 'ylim', [-2 1])
set(gca, 'xlim', [-2 2])

zeroAxes(gca);

void = semilogx(NaN, NaN, 'k.-');
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 24)
set(gca, 'ylim', [-12 -4])
set(gca, 'xlim', [0.025 32]) 
grid on
ylabel('log-likelihood')
xlabel('length scale')

clf

for i = 1:length(lengthScale)
  kern.comp{1}.inverseWidth = 1/(lengthScale(i)*lengthScale(i));
  K = kernCompute(kern, x);
  [invK, U] = pdinv(K);
  logDetK = logdet(K, U);
  ll(i) = -0.5*(logDetK + y'*invK*y + size(y, 1)*log(2*pi));
  llLogDet(i) = -.5*(logDetK+size(y, 1)*log(2*pi));
  llFit(i) = -.5*y'*invK*y;
  Kx = kernCompute(kern, x, xtest);
  ypredMean = Kx'*invK*y;
  ypredVar = kernDiagCompute(kern, xtest) - sum((Kx'*invK).*Kx', 2);

  counter = counter + 1;
  figure(counter)
  clf
  fill([xtest; xtest(end:-1:1)], ...
       [ypredMean; ypredMean(end:-1:1)] ...
         + 2*[ypredVar; -ypredVar], ...
       fillColor,'EdgeColor',fillColor)
  hold on;
  t = plot(xtest, ypredMean, 'k-');
  
  p = plot(x, y, markerType);
  set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
  set(t, 'linewidth', lineWidth);
  set(gca, 'fontname', 'times')
  set(gca, 'fontsize', 40)
  set(gca, 'ylim', [-2 1])
  
  zeroAxes(gca, 0.025, 40);
  fileName = ['demOptimiseGpTutorial' num2str(counter)];
  if exist('printDiagram') && printDiagram
    printPlot(fileName, '../tex/diagrams');
  end
  if i == length(lengthScale)
    counter = counter + 1;
    figure(counter)
    t = semilogx(lengthScale(1:i), -ll(1:i), 'k.-');
    hold on
    t = [t; semilogx(lengthScale(1:i), -llLogDet(1:i), 'k.:')];
    t = [t; semilogx(lengthScale(1:i), -llFit(1:i), 'k.--')];
    set(t, 'markersize', markerSize/2, 'lineWidth', markerWidth/2);
    set(gca, 'fontname', 'times')
    set(gca, 'fontsize', 18)
    set(gca, 'ylim', [0 15])
    set(gca, 'xlim', [0.025 32]) 
    grid on
    ylabel('negative log-likelihood')
    xlabel('length scale')
    fileName = ['demOptimiseGpTutorial' num2str(counter)];
    if exist('printDiagram') && printDiagram
      printPlot(fileName, '../tex/diagrams');
    end
  end
end
  


