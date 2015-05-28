function ll = mlpLogLikelihood(model)

% MLPLOGLIKELIHOOD Multi-layer perceptron log likelihood.
% FORMAT
% DESC computes the log likelihood of a multi-layer perceptron
% model. For single hidden layer models this is done by wrapping 
% the mlperr command. 
% ARG model : the model structure for computing the log likelihood.
% RETURN ll : the model log likelihood.
%
% SEEALSO : modelLogLikeihood, mlperr
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS


if length(model.hiddenDim) == 1
  ll = -mlperr(model, model.X, model.y);
else
  Y = mlpOut(model, model.X);
  ll = -0.5*sum(sum((model.Y - Y).^2));
end

ll = ll - size(model.X, 1)/2*log(2*pi);
