function modelRet = gpsimTest

% GPSIMTEST Test the gradients of the GPSIM model.
% FORMAT
% DESC runs some tests on the code in the GPSIM toolbox to
% test that it is working.
% RETURN model : a cell array of models used for testing.
%
% SEEALSO : modelTest
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006


% SHEFFIELDML

numGenes = 5;
numProteins = 1;
numData = 4;
times = linspace(0, 5, numData)';
y = randn(numData, numGenes);
yVar = randn(numData, numGenes);
yVar = yVar.*yVar;

numCandGenes = 3;
numCandData = 5;
timesCand = linspace(0, 5, numCandData)';
yCand = randn(numCandData, numCandGenes);
yCandVar = randn(numCandData, numCandGenes);
yCandVar = yCandVar.*yCandVar;

options = gpsimOptions;

model = gpsimCreate(numGenes, numProteins, times, y, yVar, options);

fprintf('Standard Parameter test:\n');
modelGradientCheck(model);

model = gpsimAddCandidate(model, numCandGenes, timesCand, yCand, yCandVar, options);

fprintf('Candidate Parameter test:\n');
model.type = 'gpsimCandidate';
modelGradientCheck(model);
model.type = 'gpsim';

modelRet = model;
