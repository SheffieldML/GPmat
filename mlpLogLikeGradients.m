function g = mlpLogLikeGradients(model)

% MLPLOGLIKEGRADIENTS Multi-layer perceptron gradients.
% FORMAT
% DESC computes the gradients of the log likelihood of a
% multi-layer perceptron with respect to the parameters. For
% networks with one hidden layer this is done by wrapping the
% mlpgrad command.
% ARG model : the model structure for computing the log likelihood.
% RETURN g : the gradients of the model log likelihood.
%
% SEEALSO : modelLogLikeihood, mlpgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

if length(model.hiddenDim) == 1
  g = -mlpgrad(model, model.X, model.y);
else
  error('Not yet implemented.')
end
