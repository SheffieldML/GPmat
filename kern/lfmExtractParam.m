function [param, names] = lfmExtractParam(model)

% LFMEXTRACTPARAM Extract the parameters of an LFM model.
% FORMAT
% DESC extracts the model parameters from a structure containing
% the information about a Gaussian process latent force model.
% ARG model : the model structure containing the information about
% the model.
% RETURN params : a vector of parameters from the model.
%
% SEEALSO : lfmCreate, lfmExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

if nargout>1
  [param, names] = kernExtractParam(model.kern);
else
  param = kernExtractParam(model.kern);
end

if isfield(model, 'fix')
  for i = 1:length(model.fix)
    param(model.fix(i).index) = model.fix(i).value;
  end
end
param = real(param);
