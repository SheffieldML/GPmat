function model = gpCreate(Y, X, kern, approx, numActive)

% GPCREATE Create a GP model with inducing varibles/pseudo-inputs.

% FGPLVM

if nargin < 4 
  approx = 'ftc';
end
model.m = mean(Y, 1);
model.scales = ones(1, size(Y, 2));;
Y = Y - repmat(model.m, size(Y, 1), 1);
Y = Y./repmat(model.scales, size(Y, 1), 1);
model.q = size(X, 2);
model.d = size(Y, 2);
model.N = size(Y, 1);
model.Y = Y;
model.X = X;
model.optimiser = 'conjgrad';

model.type = 'gp';
model.approx = approx;

switch approx
 case 'ftc'
  model.k = 0;
  model.X_u = [];
 case {'dtc', 'fitc', 'pitc'}
  % Sub-sample inducing variables.
  model.k = numActive;
  ind = randperm(model.N);
  ind = ind(1:model.k);
  model.X_u = model.X(ind, :);
  model.sigma2 = exp(-2);
  model.sigma2Transform = 'negLogLogit';

end
if model.k>model.N
  error('Number of active points cannot be greater than number of data.')
end
if strcmp(model.approx, 'pitc')
  numBlocks = ceil(model.N/model.k);
  startVal = 1;
  endVal = model.k;
  model.blockEnd = zeros(1, numBlocks);
  for i = 1:numBlocks
    model.blockEnd(i) = endVal;
    endVal = model.k + endVal;
    if endVal>model.N
      endVal = model.N;
    end
  end  
end
if isstruct(kern) 
  model.kern = kern;
else
  model.kern = kernCreate(model.X, kern);
end

initParams = gpExtractParam(model);
% This forces kernel computation.
model = gpExpandParam(model, initParams);

