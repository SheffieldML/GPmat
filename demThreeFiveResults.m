% DEMTHREEFIVERESULTS Plot results from the three vs five experiments.
% 
% Recreates the three versus five classification experiments
% presented in the NIPS paper.

% IVM

load demThreeFiveIvm1
clf
hold on
lineA = errorbar(probLabelled, mean(areaIVM), std(areaIVM), 'rx-');
set(lineA, 'linewidth', 2);
lineB = errorbar(probLabelled, mean(areaNCNM), std(areaNCNM), 'bo:');
set(lineB, 'linewidth', 2);
lineC = errorbar(probLabelled, mean(areaSVM), std(areaSVM), 'gx-.');
set(lineC, 'linewidth', 2);
lineD = errorbar(probLabelled, mean(areaTSVM), std(areaTSVM), 'mo--');
set(lineD, 'lineWidth', 2);
set(gca, 'xscale', 'log')
set(gca, 'xlim', [0.01 0.2])
set(gca, 'ylim', [0.75 1])
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 20)
set(gca, 'ytick', [0.7 0.8 0.9 1])
xlabel('prob. of label present')
ylabel('area under ROC curve')
