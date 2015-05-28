% DEMREGRESSIONGP Demonstrate Gaussian processes for regression.
% FORMAT
% DESC runs a simple one-D Gaussian process displaying errorbars.
%
% SEEALSO : gpCreate, demInterpolationGp
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% GP

randn('seed', 1e7)
rand('seed', 1e7)

noiseLevel = 0.2;
noiseVar = noiseLevel*noiseLevel;
% Create data set
x = linspace(-1, 1, 9)';
trueKern = kernCreate(x, 'rbf');
K = kernCompute(trueKern, x) + eye(size(x, 1))*noiseVar;
% Sample some true function values.
yTrue = gsamp(zeros(size(x))', K, 1)';

fillColor =[0.7 0.7 0.7];
markerSize = 20;
markerWidth = 6;
markerType = 'k.';
lineWidth = 2;
% Create a test set
indTrain{1} = [1 9]';
indTrain{2} = [1 5 9]';
indTrain{3} = [1 3 5 7 9]';
indTrain{4} = [1 2 3 4 5 6 7 8 9]';
figNo = 1;
for i = 0:length(indTrain)
  if i > 0
    yTrain = yTrue(indTrain{i});
    xTrain = x(indTrain{i});
    kern = kernCreate(x, 'rbf');
    % Change inverse variance (1/(lengthScale^2)))
    kern.inverseWidth = 5;
    
    xTest = linspace(-2, 2, 200)';
    
    Kx = kernCompute(kern, xTest, xTrain);
    Ktrain = kernCompute(kern, xTrain, xTrain);
    
    yPred = Kx*pdinv(Ktrain + eye(size(Ktrain))*noiseVar)*yTrain;
    yVar = kernDiagCompute(kern, xTest) - sum(Kx*pdinv(Ktrain+ eye(size(Ktrain))*noiseVar).*Kx, 2);
    ySd = sqrt(yVar);
    figure(figNo)
    clf
    fill([xTest; xTest(end:-1:1)], ...
         [yPred; yPred(end:-1:1)] ...
         + 2*[ySd; -ySd], ...
         fillColor,'EdgeColor',fillColor)
    hold on;
    h = plot(xTest, yPred, 'k-');
    set(h, 'linewidth', lineWidth)
    p = plot(xTrain, yTrain, markerType);
    set(p, 'markersize', markerSize, 'lineWidth', markerWidth);
    set(gca, 'xtick', [-2 -1 0 1 2]);
    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
    zeroAxes(gca);
    if exist('printDiagram') && printDiagram
      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
    end
    figNo = figNo + 1;
  else
    xTest = linspace(-2, 2, 200)';
    kern = kernCreate(xTest, 'rbf');
    % Change inverse variance (1/(lengthScale^2)))
    kern.inverseWidth = 5;
    
    figure(figNo)
    p = [];
    yPred = zeros(size(xTest));
    ySd = sqrt(kernDiagCompute(kern, xTest));
    fill([xTest; xTest(end:-1:1)], ...
         [yPred; yPred(end:-1:1)] ...
         + 2*[ySd; -ySd], ...
         fillColor,'EdgeColor',fillColor)
    hold on;
    h = plot(xTest, yPred, 'k-');
    
    set(h, 'linewidth', lineWidth)
    set(gca, 'xtick', [-2 -1 0 1 2]);
    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
    zeroAxes(gca);
    if exist('printDiagram') && printDiagram
      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
    end
    figNo = figNo + 1;
    
  end
  if i < length(indTrain)
    figure(figNo)
    clf
    fill([xTest; xTest(end:-1:1)], ...
         [yPred; yPred(end:-1:1)] ...
         + 2*[ySd; -ySd], ...
         fillColor,'EdgeColor',fillColor)
    hold on;
    h = plot(xTest, yPred, 'k-');
    set(h, 'linewidth', lineWidth)
    p = [p plot(x(indTrain{i+1}), yTrue(indTrain{i+1}), markerType)];
    set(p, 'markersize', markerSize, 'linewidth', markerWidth);
    set(gca, 'xtick', [-2 -1 0 1 2]);
    set(gca, 'ytick', [-3 -2 -1 0 1 2 3]);
    set(gca, 'fontname', 'times', 'fontsize', 18, 'xlim', [-2 2], 'ylim', [-3 3])
    zeroAxes(gca);
    if exist('printDiagram') && printDiagram
      printPlot(['demRegression' num2str(figNo-1)], '../tex/diagrams', '../html');
    end
    figNo = figNo + 1;
  end
end
