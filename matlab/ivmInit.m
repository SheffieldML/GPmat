function model = ivmInit(model, d)

% IVMINIT Initialise the IVM model.

% IVM

if nargin < 2
  d = model.d;
end

% Get number of data.
numData = size(model.y, 1);
numOut = size(model.y, 2);

% Initialise kernel storage.
model.kern.Kstore = zeros(numData, d);

% Initialise Indices
switch model.selectionCriterion
 case 'none'
  model.I = 1:size(model.X, 1);
 otherwise
  model.I = [];
end


% Initialise parameters
model.kern.diagK = kernDiagCompute(model.X, model.kern);

model.m = sparse(zeros(numData, numOut));
model.beta = sparse(zeros(numData, numOut));

model.varSigma = repmat(model.kern.diagK, 1, numOut);
model.mu = zeros(numData, numOut);

model.g = zeros(numData, numOut);
model.nu = zeros(numData, numOut);

% Initialise site precision and mean
switch model.noise.type  
 
 case {'probit', 'multiprobit'}

  model.u = zeros(numData, numOut);
  model.c = zeros(numData, numOut);
 
 case 'heaviside'
  model.u = zeros(numData, numOut);
  model.c = zeros(numData, numOut);

 case 'ordered'
  model.u = zeros(numData, numOut);
  model.c = zeros(numData, numOut);

  

end

model = ivmUpdateNuG(model, 1:numData);

