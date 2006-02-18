function model = gpCreate(q, d, X, y, options);

% GPCREATE Create a GP model with inducing varibles/pseudo-inputs.

% FGPLVM

if size(X, 2) ~= q
  error(['Input matrix X does not have dimension ' num2str(q)]);
end
if size(y, 2) ~= d
  error(['Input matrix y does not have dimension ' num2str(d)]);
end

if any(isnan(y)) & ~options.isMissingData
  error('NaN values in y, but no missing data declared.')
end
if options.isMissingData & options.isSpherical
  error('If there is missing data, spherical flag cannot be set.');
end

model.type = 'gp';
model.approx = options.approx;

model.learnScales = options.learnScales;
model.scaleTransform = 'negLogLogit';

model.X = X;
model.y = y;

model.q = size(X, 2);
model.d = size(y, 2);
model.N = size(y, 1);

model.optimiser = options.optimiser;

model.isMissingData = options.isMissingData;
if model.isMissingData
  for i = 1:model.d
    model.indexPresent{i} = find(~isnan(y(:, i)));
  end
end

model.isSpherical = options.isSpherical;

if ~model.isMissingData
  model.bias = mean(y);
  model.scale = ones(1, model.d);
else
  for i = 1:model.d
    if isempty(model.indexPresent{i})
      model.bias(i) = 0;
      model.scale(i) = 1;
    else
      model.bias(i) = mean(model.y(model.indexPresent{i}, i));
      model.scale(i) = 1;
    end
  end
end

model.m = gpComputeM(model);

if isstruct(options.kern) 
  model.kern = options.kern;
else
  model.kern = kernCreate(model.X, options.kern);
end



switch options.approx
 case 'ftc'
  model.k = 0;
  model.X_u = [];
 case {'dtc', 'fitc', 'pitc'}
  % Sub-sample inducing variables.
  model.k = options.numActive;
  ind = randperm(model.N);
  ind = ind(1:model.k);
  model.X_u = model.X(ind, :);
  model.beta = options.beta;
  model.betaTransform = 'negLogLogit';  
end
if model.k>model.N
  error('Number of active points cannot be greater than number of data.')
end
if strcmp(model.approx, 'pitc')
  numBlocks = ceil(model.N/model.k);
  numPerBlock = ceil(model.N/numBlocks);
  startVal = 1;
  endVal = model.k;
  model.blockEnd = zeros(1, numBlocks);
  for i = 1:numBlocks
    model.blockEnd(i) = endVal;
    endVal = numPerBlock + endVal;
    if endVal>model.N
      endVal = model.N;
    end
  end  
end

initParams = gpExtractParam(model);

% This forces kernel computation.
model = gpExpandParam(model, initParams);

