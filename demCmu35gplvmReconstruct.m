% DEMCMU35GPLVMRECONSTRUCT Reconstruct right leg and body of CMU 35.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'cmu35gplvm';

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
origBias = mean(Y);
origScale = 1./sqrt(var(Y));
%scale = ones(size(scale));
Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);

% REMOVE LEG TEST
% Indices associated with right leg.
legInd = [8:14];
% Indices associated with upper body.
bodyInd = [21:50];

for experimentNo = 1:3;

  % Load saved model.
  capName = dataSetName;
  capName(1) = upper(capName(1));
  load(['dem' capName num2str(experimentNo) '.mat']);
  startInd = 63;
  
  legErrs = fgplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale, legInd, [dataSetName ...
                      'Leg'], experimentNo);
  bodyErrs = fgplvmTaylorAngleErrors(model, Y, Ytest, startInd, origBias, origScale, bodyInd, [dataSetName ...
                      'Body'], experimentNo);

  save(['dem' capName 'Reconstruct' num2str(experimentNo) '.mat'], ...
     'legErrs', 'bodyErrs');
end
