function modelRet = lfmTest

% LFMTEST Test the gradients of the LFM model.
%
%	Description:
%
%	MODEL = LFMTEST runs some tests on the code in the LFM toolbox to
%	test that it is working.
%	 Returns:
%	  MODEL - a cell array of models used for testing.
%	
%
%	See also
%	MODELTEST


%	Copyright (c) 2007 Neil D. Lawrence



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