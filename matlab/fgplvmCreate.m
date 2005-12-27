function model = fgplvmCreate(Y, latentDim, approx, numActive, kern, ...
                                 prior, backModel, varargin)

% FGPLVMCREATE Create a GPLVM model with inducing varibles.

% FGPLVM

X = ppcaEmbed(Y, latentDim);
model = gpCreate(Y, X, kern, approx, numActive);

model.type = 'fgplvm';

if nargin < 6
elseif isstruct(prior)
  model.prior = prior;
else
  model.prior = priorCreate(prior);
  if nargin >6
    if isstruct(backModel)
      model.back = backModel;
    else
      fhandle = str2func([backModel 'Create']);
      model.back = fhandle(model.d, model.q, varargin{:});
    end
    model.back = mappingOptimise(model.back, model.Y, model.X);
  end
end

initParams = fgplvmExtractParam(model);
% This forces kernel computation.
model = fgplvmExpandParam(model, initParams);

