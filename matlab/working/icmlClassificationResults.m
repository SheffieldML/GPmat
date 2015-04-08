% ICMLCLASSIFICATIONRESULTS Plot the classificaiton results for ICML paper.

% PSIVM

fontName = 'times';
fontSize = 32;
lineWidth = 2;
markerSize = 10;
seed = 1e4;

dVals = [50 75 100 150 200 250 300];
mError = zeros(size(dVals));
stdError = zeros(size(dVals));
t = zeros(size(dVals));
for i = 1:length(dVals)
  d = dVals(i);
  load(['icmlVowel_d' num2str(d) '_seed' num2str(seed)])
  t(i) = mean(time);
  error = 100 - perCentCorrect;
  mError(i) = mean(error);
  stdError(i) = std(error);
end
linesOne = errorbar(t, mError, stdError, 'rx');
set(linesOne, 'linewidth', lineWidth);
set(linesOne, 'markersize', markerSize);
set(gca, 'xscale', 'log')

dVals = [100 200 300 400 500 700];
mError = zeros(size(dVals));
stdError = zeros(size(dVals));
t = zeros(size(dVals));
for i = 1:length(dVals)
  d = dVals(i);
  load(['icmlPsVowel_d' num2str(d) '_seed' num2str(seed)])
  t(i) = mean(time);
  error = 100 - perCentCorrect;
  mError(i) = mean(error);
  stdError(i) = std(error);
end
hold on
linesTwo = errorbar(t, mError, stdError, 'bo')
set(linesTwo, 'linewidth', lineWidth);
set(linesTwo, 'markersize', markerSize);

xlab = xlabel('time/s');
set(xlab, 'fontname', fontName, 'fontsize', fontSize)
ylab = ylabel('% error');
set(ylab, 'fontname', fontName, 'fontsize', fontSize)
set(gca,'fontname', fontName, 'fontsize', 7/8*fontSize)

set(gca, 'xlim', [50 10000])
grid on
set(gca, 'xtick', [100 1000 10000])
print -deps ../tex/diagrams/classificationComparison.eps