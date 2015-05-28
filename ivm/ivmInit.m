function model = ivmInit(model, d)

% IVMINIT Initialise the IVM model.
% FORMAT
% DESC sets up some initial matrices and vectors in the IVM
% representation in preparation for learning.
% ARG model : the model to be set up.
% ARG d : optional value of the active set size (if not given,
% model.d is used).
% RETURN model : model with various matrices set up for learning.
%
% SEEALSO : ivmOptimiseIvm, ivmSelectPoints
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

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

model.terminate = 0;

% Initialise parameters
model.kern.diagK = kernDiagCompute(model.kern, model.X);

model.m = sparse(zeros(numData, numOut));
model.beta = sparse(zeros(numData, numOut));

model.varSigma = repmat(model.kern.diagK, 1, numOut);
model.mu = zeros(numData, numOut);

model.g = zeros(numData, numOut);
model.nu = zeros(numData, numOut);

%/~
% % Initialise site precision and mean
% switch model.noise.type  
 
%  case {'probit', 'multiprobit'}

%   model.u = zeros(numData, numOut);
%   model.c = zeros(numData, numOut);
 
%  case 'heaviside'
%   model.u = zeros(numData, numOut);
%   model.c = zeros(numData, numOut);

%  case 'ordered'
%   model.u = zeros(numData, numOut);
%   model.c = zeros(numData, numOut);

  

% end
%~/
model = ivmUpdateNuG(model, 1:numData);

