function model = dnetCreate(inputDim, outputDim, Y, options)

% DNETCREATE Density network model.
%
%	Description:
%
%	MODEL = DNETCREATE(INPUTDIMENSION, OUTPUTDIM, Y, OPTIONS) creates a
%	structure for a density network.
%	 Returns:
%	  MODEL - model structure containing the density network specified.
%	 Arguments:
%	  INPUTDIMENSION - dimension of latent data.
%	  OUTPUTDIM - dimension of observed data.
%	  Y - the data to be modelled in design matrix format (as many rows
%	   as there are data points).
%	  OPTIONS - options structure as returned by dnetCreate.
%	
%
%	See also
%	DNETOPTIONS, MLPCREATE, RBFCREATE, KBRCREATE


%	Copyright (c) 2008 Neil D. Lawrence



model.type = 'dnet';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end

if isstr(options.initX)
  initFunc = str2func([options.initX 'Embed']);
  model.X = initFunc(Y, inputDim);
else
  if size(options.initX, 1) == size(Y, 1) ...
        & size(options.initX, 2) == inputDim
    model.X = options.initX;
  else
    error('options.initX not in recognisable form.');
  end
end


model.q = inputDim;
model.d = outputDim;
model.N = size(Y, 1);

model.mapping = modelCreate(options.mappingType, inputDim, ...
                            outputDim, options.mappingOptions);
model.M = options.M;
model.y = Y;
model.w = repmat(1/model.M, model.N, model.M);
model.betaTransform =  optimiDefaultConstraint('positive');  

model.beta = exp(-2);
model.grid = options.grid;

if ~isempty(model.grid)
  minVals = 1.2*min(model.X);
  maxVals = 1.2*max(model.X);
  x = linspace(minVals(1), maxVals(1), model.grid(1));
  y = linspace(minVals(2), maxVals(2), model.grid(2));
  [mx, my] = meshgrid(x, y);
  for i = 1:size(model.grid)
    model.X_u = [mx(:) my(:)];
  end
else
  model.X_u = randn(model.M, model.q);
end

logw = -dist2(model.X, model.X_u);
logw = logw -repmat(max(logw, [], 2), 1, model.M);
model.w = exp(logw)./repmat(sum(exp(logw), 2), 1, model.M);

model.basisStored = options.basisStored;

model.alpha = options.alpha;

params = dnetExtractParam(model);
model = dnetExpandParam(model, params);

if model.basisStored
  model = dnetUpdateOutputWeights(model);
  model = dnetUpdateBeta(model);  
end

model.numParams = model.mapping.numParams + 1;
model.outputDim = model.d;
model.inputDim = model.q;