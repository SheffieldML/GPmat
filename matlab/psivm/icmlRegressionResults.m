% ICMLREGRESSIONRESULTS Prepare comparision plots of the IVM vs sub-sampling.

% PSIVM

markerSize = 10;
fontName = 'times';
fontSize = 32;
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
for i = 4:length(dVec)
  llMeanIVM(i) = mean(llIVM(i, :));
  llVarIVM(i) = var(llIVM(i, :));
  timeMeanIVM(i) = mean(timeIVM(i, :));
end

llModels = psivm(testX, testY, 'rbf', 'gaussian', 'none', []);
llModels.lntheta = log(thetaConstrain([1 1 100 0]));
llModels = psivmInit(llModels);
bestLl = -pskernelObjective(log(thetaConstrain([1 1 100 0])), llModels, ...
			    0); 
numTestTasks = length(testX);
colordef white

linesTime = errorbar([timeMeanSub timeMeanIVM], (bestLl-[llMeanSub llMeanIVM])/numTestTasks, ...
		    sqrt([llVarSub llVarIVM]/numTestTasks^2));
set(linesTime, 'linewidth', 2);
set(linesTime(end/2+1:3*end/4), ...
    'linestyle', 'none', ...
    'marker', 'x', ...
    'markersize', markerSize)
set(linesTime(3*end/4+1:end), ...
    'linestyle', 'none', ...
    'marker', 'o', ...
    'markersize', markerSize)

xlab = xlabel('time/s');
set(xlab, 'fontname', fontName,  'fontsize', fontSize)
ylab = ylabel('KL');
set(ylab, 'fontname', fontName, 'fontsize', fontSize)
grid on
set(gca,'fontname', fontName, 'fontsize', fontSize*7/8)
set(gca, 'ylim', [-1 20])
set(gca, 'xlim', [0 800])
print -deps ../tex/diagrams/fullPlot.eps

figure
linesSamps = errorbar([sampsVec' dVec'], (bestLl-[llMeanSub llMeanIVM])/numTestTasks, ...
		    sqrt([llVarSub llVarIVM]/numTestTasks^2));
set(gca, 'ylim', [-1 20])
set(gca, 'xlim', [0 800])

%set(linesSamps(1), 'linestyle', 'none', 'marker', 'x', 'markersize', ...
%		  markerSize)
%set(linesSamps(2), 'linestyle', 'none', 'marker', 'o', 'markersize', ...
%		  markerSize)
set(linesSamps, 'linewidth', 2);
set(linesSamps, 'linewidth', 2);
set(linesSamps(end/2+1:3*end/4), ...
    'linestyle', 'none', ...
    'marker', 'x', ...
    'markersize', markerSize)
set(linesSamps(3*end/4+1:end), ...
    'linestyle', 'none', ...
    'marker', 'o', ...
    'markersize', markerSize)

xlab = xlabel('time/s');
set(xlab, 'fontname', fontName,  'fontsize', fontSize)
ylab = ylabel('KL');
set(ylab, 'fontname', fontName, 'fontsize', fontSize)
grid on
set(gca,'fontname', fontName, 'fontsize', fontSize*7/8)
set(gca, 'ylim', [-1 20])
%set(gca, 'xlim', [0 800])
