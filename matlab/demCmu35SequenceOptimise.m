% DEMCMU35SEQUENCEOPTIMISE 

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'cmu35gplvm';

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
skel = acclaimReadSkel('35.asf');
[tmpchan, skel] = acclaimLoadChannels('35_01.amc', skel);


% REMOVE LEG TEST
% Indices associated with right leg.
legInd = [8:14];
% Indices associated with upper body.
bodyInd = [21:50];

% Load saved model.
capName = dataSetName;
capName(1) = upper(capName(1));

for experimentNo = 1:3;
  load(['dem' capName num2str(experimentNo) '.mat']);
  for missing = 1:2
    if missing==1
      type = 'Leg';
      missingInd = legInd;
    else
      type = 'Body';
      missingInd = bodyInd;
    end
    if exist(['dem' capName 'Yvals' type num2str(experimentNo) '.mat'], 'file')==2
      load(['dem' capName 'Yvals' type num2str(experimentNo) '.mat'])
    else
      disp('Reconstructing ...');
      startInd = 63;
      
      Ytest(startInd:end, missingInd) = NaN;
      model = gpComputeAlpha(model);
      ll = [];
      for j = 1:size(model.X_u, 1)
        for i = 1:size(Ytest, 1);
          ll(i, j) = fgplvmPointLogLikelihood(model, model.X_u(j, :), ...
                                              Ytest(i, :));
        end
      end
      [void, ind] = max(ll, [], 2);
      Xinit = model.X_u(ind, :);
      Xpred = fgplvmOptimiseSequence(model, Xinit, Ytest, 1, 1000);
      Xpred = Xpred(startInd:end, :);
      Ytest = Ytest(startInd:end, :);
      Ypred = gpOut(model, Xpred);
      
      
      save(['dem' capName 'Yvals' type num2str(experimentNo) '.mat'], ...
           'Xpred', 'Ypred');
    end
  end
end
