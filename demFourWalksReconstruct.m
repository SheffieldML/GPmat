% DEMFOURWALKSRECONSTRUCT Reconstruct right leg of CMU 35.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'fourWalks';
experimentNo = 1;

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
YtrueTest = Ytest;
indKeepDim = find(var(Y)>1e-6);


% Load saved model.
capName = dataSetName;;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat']);
  model = gpComputeAlpha(model);

startObserved = [1:4];
endObserved = [36:50];
Xinit = zeros(size(Ytest, 1), 3);
Xtest = zeros(size(Ytest, 1), 3);
for k = 1
%for k = 1:8
  st = 50*(k-1)+1;
  en1 = 50*(k-1)+36;
  en2 = 50*k;
  startObserved = [st:st+3];
  endObserved = [en1:en2];
  observedIndex = [startObserved endObserved];
  missingIndex = (max(startObserved)+1):(min(endObserved)-1);
  Ytest(missingIndex, :) = NaN;

  if k == 1
    [Xinit(st:en2, :), llo] = fgplvmViterbiSequence(model, model.X, ...
                                                    Ytest(st:en2, :));
  else
    Xinit(st:en2, :) = fgplvmViterbiSequence(model, Ytest(st:en2, :), llo);
  end

  %   % Use forward prediction dynamics to initalise
  %   for i = min(missingIndex):max(missingIndex);
  %     Xinit(i, :) = gpPosteriorMeanVar(model.dynamics, Xinit(i-1, :));
  %     if model.dynamics.diff
  %       Xinit(i, :) = Xinit(i-1, :) + Xinit(i, :);
  %     end
  %   end

  Xtest(st:en2, :) = fgplvmOptimiseSequence(model, Xinit(st:en2, :), ...
                                            Ytest(st:en2, :), 1, 5000);
  Ypred(st:en2, :) = gpOut(model, Xtest(st:en2, :));

  diff = Ypred(missingIndex, indKeepDim)- YtrueTest(missingIndex, indKeepDim);
  diff = diff(:);
  diff = diff.*diff;
  rmse(k) = sqrt(mean(diff(:)));
end

%save demFourWalksReconstruct Xinit Xtest Ypred
