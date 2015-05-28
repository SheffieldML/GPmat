function g = dnetOutputGradX(model, X)

% DNETOUTPUTGRADX Evaluate derivatives of DNET model outputs with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an DNET model with
% respect to the inputs to the model. 
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs.
%
% SEEALSO : dnetOutputGrad, modelOutputGradX
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

g = modelOutputGradX(model.mapping, X);
