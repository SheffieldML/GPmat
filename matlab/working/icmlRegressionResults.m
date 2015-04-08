% ICMLREGRESSIONRESULTS Prepare comparision plots of the IVM vs sub-sampling.

% PSIVM

load psRegressionData.mat 
load icmlPsRegressionResults.mat

llMeanSub = repmat(NaN, length(sampsVec), 1);
llVarSub = repmat(NaN, length(sampsVec), 1);
timeMeanSub = repmat(NaN, length(sampsVec), 1);
for i = 1:length(sampsVec)
    llMeanSub(i) = mean(llSub(i, :));
    llVarSub(i) = var(llSub(i, :));
    timeMeanSub(i) = mean(timeSub(i, :));
end
llMeanIVM = repmat(NaN, length(dVec), 1);
llVarIVM = repmat(NaN, length(dVec), 1);
timeMeanIVM = repmat(NaN, length(dVec), 1);
% Ignore the first two d values as they have v big errorbars and confuse
% the plot
for i = 3:length(dVec)
    llMeanIVM(i) = mean(llIVM(i, :));
    llVarIVM(i) = var(llIVM(i, :));
    timeMeanIVM(i) = mean(timeIVM(i, :));
end
llModels = psivm(testX, testY, 'rbf', 'gaussian', 'none', []);
llModels.lntheta = log(thetaConstrain([1 1 100 0]));
llModels = psivmInit(llModels);
bestLl = -pskernelObjective(log(thetaConstrain([1 1 100 0])), llModels, 0); 
numTestTasks = length(testX);
colordef white
clf
%plot(timeMeanSub, bestLl-llMeanSub, 'rx--')
%hold on
linesSub = errorbar([timeMeanSub timeMeanIVM], (bestLl-[llMeanSub llMeanIVM])/numTestTasks, ...
		    sqrt([llVarSub llVarIVM]/numTestTasks^2));
%set(gca, 'ylim', [-1 5])
%set(gca, 'xlim', [0 1750])
set(linesSub, 'linewidth', 2);
markerSize = 10;
set(linesSub(end/2+1:3*end/4), 'linestyle', 'none', 'marker', 'x', 'markersize', ...
		  markerSize)
set(linesSub(3*end/4+1:end), 'linestyle', 'none', 'marker', 'o', 'markersize', ...
		  markerSize)

xlab = xlabel('time/s');
set(xlab, 'fontname', 'times', 'fontsize', 32)
ylab = ylabel('KL');
set(ylab, 'fontname', 'times', 'fontsize', 32)
%hold on
%linesIVM = errorbar(timeMeanIVM, bestLl -llMeanIVM, sqrt(llVarIVM), 'bo-')
grid on
%set(linesIVM, 'linewidth', 2);
set(gca,'fontname', 'times', 'fontsize', 28)
set(gca, 'ylim', [-1 20])
set(gca, 'xlim', [0 800])
%set(gca, 'xscale', 'log')
print -deps ../tex/diagrams/fullPlot.eps
%set(gca, 'xlim', [50 400])
%set(gca, 'ylim', [-1 5])
%print -deps ../tex/diagrams/zoomPlot.eps
