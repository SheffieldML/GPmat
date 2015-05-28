function modelRet = lfmTest

% LFMTEST Test the gradients of the LFM model.
% FORMAT
% DESC runs some tests on the code in the LFM toolbox to
% test that it is working.
% RETURN model : a cell array of models used for testing.
%
% SEEALSO : modelTest
%
% COPYRIGHT : Neil D. Lawrence, 2007


% KERN

numDisplacements = 2;
numForces = 1;
numData = 10;
times = linspace(0, 5, numData)';
y = randn(numData, numDisplacements);

options = lfmOptions;

model = lfmCreate(numDisplacements, numForces, times, y, options);

fprintf('Standard Parameter test:\n');
modelGradientCheck(model);

modelRet = model;
