function g = linearOutputGradX(model, X)

% LINEAROUTPUTGRADX Evaluate derivatives of linear model outputs with respect to inputs.
% FORMAT
% DESC returns the derivatives of the outputs of an LINEAR model with
% respect to the inputs to the model. 
% ARG model : the model for which the derivatives will be computed.
% ARG X : the locations at which the derivatives will be computed.
% RETURN g : the gradient of the output with respect to the inputs, in
%
% SEEALSO : linearOutputGrad, modelOutputGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

g = repmat(shiftdim(model.W, -1), [size(X, 1) 1 1]);
