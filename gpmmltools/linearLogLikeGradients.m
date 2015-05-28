function g = linearLogLikeGradients(model)

% LINEARLOGLIKEGRADIENTS Linear model gradients.
% FORMAT
% DESC computes the gradients of the log likelihood of a
% linear model with respect to the parameters.
% ARG model : the model structure for computing the log likelihood.
% RETURN g : the gradients of the model log likelihood.
%
% SEEALSO : modelLogLikeihood, lineargrad
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

Xo = [model.X ones(size(model.X, 1), 1)];
W = [model.W; model.b];
G = -model.beta*(Xo'*Xo*W - Xo'*model.y);
gW = G(1:end-1, :);
gb = G(end, :);
g = [gW(:)' gb];
