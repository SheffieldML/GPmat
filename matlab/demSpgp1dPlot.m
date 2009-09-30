% DEMSPGP1DPLOT Plot results from 1-D sparse GP.

% GP

fillColor = [0.7 0.7 0.7];
xTest = linspace(-1.5, 1.5, 200)';
[mu, varSigma] = gpPosteriorMeanVar(model, xTest);
%[mu, varSigma] = gpVarPosteriorMeanVar(model, xTest);


figure
fill([xTest; xTest(end:-1:1)], ...
     [mu; mu(end:-1:1)] ...
     + 2*[sqrt(varSigma); -sqrt(varSigma(end:-1:1))], ...
     fillColor,'EdgeColor',fillColor)
hold on;
plot(X, y, 'k.');
a = plot(xTest, mu, 'k-');
%/~
%a = [a plot(xTest, mu+2*sqrt(varSigma), 'b--')];

%a = [a plot(xTest, mu-2*sqrt(varSigma), 'b--')];
%~/
if isfield(model, 'X_u') && ~isempty(model.X_u)
  b = plot(model.X_u, -ones(size(model.X_u)), 'bx');
  set(b, 'linewidth', 2)
  set(b, 'markersize', 10);
end
set(gca, 'ylim', [-1 2])
set(gca, 'xlim', [-1.5 1.5])
set(a, 'linewidth', 2);
zeroAxes(gca, [], 10, 'arial')
if exist('printDiagram') && printDiagram
  printPlot(fileName, '../tex/diagrams', '../html');
end
