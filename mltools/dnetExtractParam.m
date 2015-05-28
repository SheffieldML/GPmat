function [params, names] = dnetExtractParam(model)

% DNETEXTRACTPARAM Extract weights and biases from an DNET.
% FORMAT
% DESC returns a vector of all the parameters from a density network.
% ARG model : the model from which we wish to extract the weights
% and biases.
% RETURN params : vector of all the weights and biases returned by
% the model. The structure is governed by dnetpak.
% RETURN names : optional additional returned cell array of the
% names of the parameters.
%
% SEEALSO : dnetCreate, dnetExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

if nargout > 1
  [params, names] = modelExtractParam(model.mapping);
  names{end+1} = 'Output Precision';
else
  params = modelExtractParam(model.mapping);
end
func = str2func([model.betaTransform 'Transform']);
params(end+1) = func(model.beta, 'xtoa');
