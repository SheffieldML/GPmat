function model = fgplvmCreate(q, d, Y, options)

% FGPLVMCREATE Create a GPLVM model with inducing variables.
% FORMAT
% DESC creates a GP-LVM model with the possibility of using
% inducing variables to speed up computation.
% ARG q : dimensionality of latent space.
% ARG d : dimensionality of data space.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned from
% FGPLVMOPTIONS. This structure determines the type of
% approximations to be used (if any).
% RETURN model : the GP-LVM model.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2010
%
% SEEALSO : modelCreate, fgplvmOptions

% FGPLVM

if size(Y, 2) ~= d
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end

if isstr(options.initX)
  initFunc = str2func([options.initX 'Embed']);
  X = initFunc(Y, q);
else
  if size(options.initX, 1) == size(Y, 1) ...
        & size(options.initX, 2) == q
    X = options.initX;
  else
    error('options.initX not in recognisable form.');
  end
end


model = gpCreate(q, d, X, Y, options);

model.type = 'fgplvm';

if isstruct(options.prior)
  model.prior = options.prior;
else
  if ~isempty(options.prior)
    model.prior = priorCreate(options.prior);
  end
end

if isstruct(options.inducingPrior)
  model.inducingPrior = options.inducingPrior;
else
  if ~isempty(options.inducingPrior)
    model.inducingPrior = priorCreate(options.inducingPrior);
  end
end

if isfield(options, 'back') & ~isempty(options.back)
  if isstruct(options.back)
    model.back = options.back;
  else
    if ~isempty(options.back)
      model.back = modelCreate(options.back, model.d, model.q,options.backOptions);
      if(isfield(options.backOptions,'indexOut')&&~isempty(options.backOptions.indexOut))
        model.back.indexOut = options.backOptions.indexOut;
      else
        model.back.indexOut = 1:1:model.q;
      end
    end
  end
  if options.optimiseInitBack
    % Match back model to initialisation.
    model.back = mappingOptimise(model.back, model.y, model.X(:,model.back.indexOut));
  end
  % Now update latent positions with the back constraints output.
  model.X(:,model.back.indexOut) = modelOut(model.back, model.y);
else
  model.back = [];
end

model.constraints = {};

model.dynamics = [];

initParams = fgplvmExtractParam(model);
model.numParams = length(initParams);
% This forces kernel computation.
model = fgplvmExpandParam(model, initParams);

