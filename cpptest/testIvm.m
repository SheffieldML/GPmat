% probit

% Gaussian

% Sample a regression data-set.
[X, y] = ivmLoadData('regressionOne');
noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbf', 'lin', 'bias', 'white'};

options = ivmOptions;
dVal = 50;
% Initialise the IVM model.

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model = ivmOptimiseIvm(model, 0);
[kernInit, noiseInit, ivmInfoInit]  = ivmDeconstruct(model);


save testGaussian X y ivmInfoInit kernInit noiseInit

% Probit

% Sample a classification data-set.
[X, y] = ivmLoadData('classificationTwo');
noiseModel = 'probit';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbf', 'lin', 'bias', 'white'};

options = ivmOptions;
dVal = 50;
% Initialise the IVM model.

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model = ivmOptimiseIvm(model, 0);
[kernInit, noiseInit, ivmInfoInit]  = ivmDeconstruct(model);


save testProbit X y ivmInfoInit kernInit noiseInit

% Ncnm

% Sample a classification data-set.
[X, y] = ivmLoadData('classificationTwo');
missing = rand(size(y))>.1;
y(missing)=NaN;
noiseModel = 'ncnm';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbf', 'lin', 'bias', 'white'};

options = ivmOptions;
dVal = 50;
% Initialise the IVM model.

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model = ivmOptimiseIvm(model, 0);
[kernInit, noiseInit, ivmInfoInit]  = ivmDeconstruct(model);


save testNcnm X y ivmInfoInit kernInit noiseInit
