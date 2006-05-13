function ll = mlpLogLikelihood(model)

% MLPLOGLIKELIHOOD Multi-layer perceptron log likelihood.
% FORMAT
% DESC computes the log likelihood of a multi-layer perceptron
% model. This is done by wrapping the mlperr command. 
% ARG model : the model structure for computing the log likelihood.
% RETURN ll : the model log likelihood.
%
% SEEALSO : modelLogLikeihood, mlperr
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

ll = -mlperr(model, model.X, model.y);
