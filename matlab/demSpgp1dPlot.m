% DEMSPGP1DPLOT Plot results from 1-D sparse GP.

% GP

xTest = linspace(-1.5, 1.5, 200)';
[mu, varSigma] = gpPosteriorMeanVar(model, xTest);

figure
plot(X, y, 'r.');
hold on
a = plot(xTest, mu, 'b-');
a = [a plot(xTest, mu+2*sqrt(varSigma), 'b--')];

a = [a plot(xTest, mu-2*sqrt(varSigma), 'b--')];
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
  fileName = ['dem' capName num2str(experimentNo)];
  print('-depsc', ['../tex/diagrams/' fileName])
  set(a, 'linewidth', 1);

  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
%  fontsize = get(gca, 'fontsize');
%  set(gca, 'fontsize', fontsize/2);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
  set(gcf, 'paperposition', origpos)
end