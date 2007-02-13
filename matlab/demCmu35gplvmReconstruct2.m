% DEMCMU35GPLVMRECONSTRUCT Reconstruct right leg of CMU 35.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'cmu35gplvm';
experimentNo = 2;

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
origBias = mean(Y);
origScale = 1./sqrt(var(Y));
%scale = ones(size(scale));
Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);
origYtest = Ytest;
startInd = 63;
Ytest = Ytest(startInd:end, :);
% Indices associated with right leg.
legInd = [8:14];
leglessInd = [1:7 15:size(Ytest, 2)];
YtrainLeg = Y(:, leglessInd);
YtestLeg = Ytest(:, leglessInd);

dists = dist2(YtrainLeg, YtestLeg);
[void, bestIndLeg] = min(dists);
lenVal = size(Ytest, 1);
nNerrLeg = (Ytest(:, legInd)-Y(bestIndLeg, legInd))./repmat(origScale(legInd), lenVal, 1);
nNerrLeg = nNerrLeg*180/pi;
nNerrLeg = mean(mean((nNerrLeg.*nNerrLeg)));


% Load saved model.
capName = dataSetName;;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat']);

YtrueTest = Ytest;
Ytest(:, legInd) = NaN;
model = gpComputeAlpha(model);
for j = 1:size(model.X_u, 1)
  for i = 1:size(Ytest, 1);
    ll(i, j) = fgplvmPointLogLikelihood(model, model.X_u(j, :), ...
                                        Ytest(i, :));
  end
end
[void, ind] = max(ll, [], 2);
Xinit = model.X_u(ind, :);
Xtest = fgplvmOptimiseSequence(model, Xinit, Ytest, 1, 1000);
Ypred = gpOut(model, Xtest);
% De-normalise the predictions.
errLeg = (YtrueTest(:, legInd) - Ypred(:, legInd))./repmat(origScale(1, legInd), size(YtrueTest, 1), 1);
% Convert to degrees.
errLeg = errLeg*180/pi;
% Compute average of mean square error.
errLeg = mean(mean((errLeg.*errLeg)));
save(['dem' capName 'Reconstruct.mat'], 'Xtest', 'Ytest')

errLeg = (YtrueTest(:, legInd) - Ypred(:, legInd))...
         ./repmat(origScale(1, legInd), size(YtrueTest, 1), 1);
errLeg = errLeg*180/pi;
load cmu35TaylorScaleBias
errLeg = errLeg.*repmat(scale(1, legInd+9), size(YtrueTest, 1), 1);
errLeg = sum(sum(errLeg.*errLeg))/7;

nNerrLeg = (YtrueTest(:, legInd) - Y(bestIndLeg, legInd))...
         ./repmat(origScale(1, legInd), size(YtrueTest, 1), 1);
nNerrLeg = nNerrLeg*180/pi;
load cmu35TaylorScaleBias
nNerrLeg = nNerrLeg.*repmat(scale(1, legInd+9), size(YtrueTest, 1), 1);
nNerrLeg = sum(sum(nNerrLeg.*nNerrLeg))/7;

plotRange = 8:14;
for plotNo = plotRange
  figNo = plotNo - min(plotRange) + 1;
  figure(figNo)
  clf
  lin = plot(1:size(origYtest, 1), origYtest(:, plotNo), '-')
  hold on
  lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Y(bestIndLeg, plotNo)], ':')];
  lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
  ax = gca;
  xlabel('Frame')
  ylabel('Normalized joint angle')
  set(ax, 'fontname', 'arial');
  set(ax, 'fontsize', 20);
  set(lin, 'lineWidth', 2);
  fileName = ['dem' capName 'Reconstruct' num2str(figNo)];
  print('-depsc', ['../tex/diagrams/' fileName])
  print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
  
  % make smaller for PNG plot.
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  fontsize = get(gca, 'fontsize');
  set(gca, 'fontsize', fontsize/2);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
  set(gcf, 'paperposition', origpos);
end
