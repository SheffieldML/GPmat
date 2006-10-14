function model = gpExpandParam(model, params)

% GPEXPANDPARAM Expand a parameter vector into a GP model.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG params : a vector of parameters for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameters.
% 
% SEEALSO : gpCreate, gpExtractParam, modelExtractParam, gpUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM


if strcmp(model.approx, 'ftc') | model.fixInducing
  endVal = 0;
else
  startVal = 1;
  endVal = model.k*model.q;
  model.X_u = reshape(params(startVal:endVal), model.k, model.q);
end
startVal = endVal +1;
endVal = endVal + model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

if model.learnScales
  startVal = endVal + 1;
  endVal = endVal + model.d;
  fhandle = str2func([model.scaleTransform 'Transform']);
  model.scale = fhandle(params(startVal:endVal), 'atox');
  model.m = gpComputeM(model);
end

% Record the total number of parameters.
model.nParams = endVal;

% Update the kernel representations.
switch model.approx
 case 'ftc'
  model = gpUpdateKernels(model, model.X, model.X_u);
 case {'dtc', 'fitc', 'pitc'}
  model = gpUpdateKernels(model, model.X, model.X_u, params(end));
 otherwise
  error('Unknown approximation type.')
end

% Update the vector 'alpha' for computing posterior mean.
if isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end



