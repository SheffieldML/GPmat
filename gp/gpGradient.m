function g = gpGradient(params, model)

% GPGRADIENT Gradient wrapper for a GP model.
% FORMAT 
% DESC wraps the log likelihood gradient function to return the
% gradient of the negative of the log likelihood. This can then be
% used in, for example, NETLAB, minimisation tools.
% ARG params : the parameters of the model.
% ARG model : the model for which gradients will be computed.
% RETURN g : the returned gradient of the negative log likelihood
% for the given parameters.
%
% SEEALSO : scg, conjgrad, gpCreate, gpGradient, gpLogLikeGradient, gpOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% GP

model = gpExpandParam(model, params);
g = - gpLogLikeGradients(model);
