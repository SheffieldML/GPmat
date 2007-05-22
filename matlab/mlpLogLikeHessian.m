function g = mlpLogLikeHessian(model)

% MLPLOGLIKEHESSIAN Multi-layer perceptron Hessian.
% FORMAT
% DESC computes the Hessian of the log likelihood of a
% multi-layer perceptron with respect to the parameters. For
% networks with a single hidden layer this is
% done by wrapping the mlpgrad command.
% ARG model : the model structure for computing the log likelihood.
% RETURN g : the Hessian of the model log likelihood.
%
% SEEALSO : modelLogLikeihood, mlpgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

if length(model.hiddenDim)==1
  g = -mlphess(model, model.X, model.y);
else
  error('Hessian not yet available for this model.')
end
