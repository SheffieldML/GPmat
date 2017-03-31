function [errStruct] = fgplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale,...
                                              missingInd, name, expNo);

% FGPLVMTAYLORANGLEERRORS Helper function for computing angle errors for CMU 35 data.

% FGPLVM

origYtest = Ytest;

YtrueTest = Ytest;

capName = name;
capName(1) = upper(capName(1));

temp = ones(1, size(Ytest, 2));
temp(missingInd) = 0;
presentInd = find(temp);
YtrainNn = Y(:, presentInd);
YtestNn = origYtest(startInd:end, presentInd);

% First nearest neighbour
dists = dist2(YtrainNn./(repmat(origScale(presentInd), size(YtrainNn, 1), 1)), ...
              YtestNn./(repmat(origScale(presentInd), size(YtestNn, 1), 1)));
[void, bestIndNn] = min(dists);
lenVal = size(Ytest, 1);
err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
         ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
err = err*180/pi;
errStruct.angleErrorNnScaled = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorNnScaled = sum(sum(err.*err))/length(missingInd);

% Second nearest neighbour
dists = dist2(YtrainNn, ...
              YtestNn);
[void, bestIndNn] = min(dists);
lenVal = size(Ytest, 1);
err = (YtrueTest(startInd:end, missingInd) - Y(bestIndNn, missingInd))...
         ./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
err = err*180/pi;
errStruct.angleErrorNn = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorNn = sum(sum(err.*err))/length(missingInd);



Ytest(startInd:end, missingInd) = NaN;
model = gpComputeAlpha(model);
for j = 1:size(model.X_u, 1)
  for i = 1:size(Ytest, 1);
    ll(i, j) = fgplvmPointLogLikelihood(model, model.X_u(j, :), ...
                                        Ytest(i, :));
  end
end
[void, ind] = max(ll, [], 2);
Xinit = model.X_u(ind, :);
try
  Xpred = fgplvmOptimiseSequence(model, Xinit, Ytest, 1, 1000);
catch e
  % Warning: This is a work-around because fgplvmOptimiseSequence seems to be missing!! The following doesn not take sequential correlations into account!
  warning('fgplvmOptimiseSequence cannot be used... reverting to fgplvmOptimisePoint.')
  Xpred = fgplvmOptimisePoint(model, Xinit, Ytest, 1, 1000);
end

Xpred = Xpred(startInd:end, :);
Ytest = Ytest(startInd:end, :);
Ypred = gpOut(model, Xpred);
% De-normalise the predictions.
err = (YtrueTest(startInd:end, missingInd) - Ypred(:, ...
                                                  missingInd))./repmat(origScale(1, missingInd), size(YtrueTest, 1)-startInd+1, 1);
% Convert to degrees.
err = err*180/pi;
% Compute average of mean square error.
errStruct.angleErrorGplvm = sqrt(mean(mean((err.*err))));
load cmu35TaylorScaleBias
err = err.*repmat(scale(1, missingInd+9), size(YtrueTest, 1)-startInd+1, 1);
errStruct.taylorErrorGplvm = sum(sum(err.*err))/length(missingInd);


plotRange = missingInd;
colordef white
for plotNo = plotRange
  figNo = plotNo - min(plotRange) + 1;
  figure(figNo)
  clf
  lin = plot(1:size(origYtest, 1), origYtest(:, plotNo), '-')
  hold on
  lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Y(bestIndNn, plotNo)], ':')];
  lin = [lin; plot(1:size(origYtest, 1), [origYtest(1:startInd-1, plotNo); Ypred(:, plotNo)], '--')];
  ax = gca;
  xlabel('Frame')
  ylabel('Normalized joint angle')
  set(ax, 'fontname', 'arial');
  set(ax, 'fontsize', 20);
  set(lin, 'lineWidth', 2);
  fileName = ['dem' capName 'Reconstruct' num2str(expNo) '_' num2str(plotNo)];
%  print('-depsc', ['../tex/diagrams/' fileName])
%  print('-deps', ['../tex/diagrams/' fileName 'NoColour'])
  %plot2svg(['../html/' fileName '.svg'])
  
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
  %print('-dpng', ['../html/' fileName])
  %plot2svg(['../html/' fileName '.svg'])
  set(gcf, 'paperposition', origpos);
end
