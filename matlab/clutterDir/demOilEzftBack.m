% DEMOILEZFTBACK Model the oil data with the EZFT algorithm and back constraints.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 1;

% load data
[Y, lbls] = gplvmLoadData(dataSetName);


% Set all points as active
numActive = 100;

% Optimise kernel parameters and active set jointly.

% Initialise X with PCA.
X = gplvmPcaInit(Y, 2);

% Plot the intialisation.
symbols = getSymbols(3);
figure, hold on
for i = 1:size(X, 1)
  labelNo = find(lbls(i, :));
  plot(X(i, 1), X(i, 2), symbols{labelNo})
end
model.type = 'fgplvm';
model.q = size(X, 2);
model.d = size(Y, 2);
numHidden = 20;
model.N = size(Y, 1);
model.k = numActive;
model.sigma2 = 0.001;
model.sigma2Transform = 'negLogLogit';
kernType = {'rbf', 'bias', 'white'};

model.Y = Y;
model.X = X;
ind = randperm(model.N);
ind = ind(1:numActive);
model.X_I = X(ind, :);

% Back constrained model.
%model.back = mlpCreate(model.d, model.q, numHidden, 'linear');
%model.back = linearCreate(model.d, model.q);
model.prior.type = 'gaussian';
model.prior = priorParamInit(model.prior);
model.prior.precision = 1;

model.kern = kernCreate(model.X_I, kernType);

initParams = gplvmEzftExtractParam(model);
% This forces kernel computation.
model = gplvmEzftExpandParam(model, initParams);

%if nargin < 4
  iters = 100;
%  if nargin < 3
    display = 1;
%    if nargin < 2
      prior = [];
%    end
%  end
%end

options = zeros(1, 18);
if display
  options(1) = 1;
  if length(initParams) <= 20
    options(9) = 1;
  end
end
options(14) = iters;

params = scg('gplvmEzftObjective', initParams,  options, ...
		    'gplvmEzftGradient', model);

model = gplvmEzftExpandParam(model, params);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% Load the results and display dynamically.
%gplvmResultsDynamic(dataSetName, experimentNo, 'vector')

% Load the results and display statically.
% gplvmResultsStatic(dataSetName, experimentNo, 'vector')

% Load the results and display as scatter plot
% gplvmResultsStatic(dataSetName, experimentNo, 'none')
