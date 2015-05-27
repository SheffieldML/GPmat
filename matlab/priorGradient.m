function g = priorGradient(prior, params)

% PRIORGRADIENT Gradient of the prior with respect to its variables
% FORMAT
% DESC wrapper function that computes the gradient of the prior with 
% respect to its parameters.
% ARG prior : the prior structure for which the gradients are being
% computed.
% ARG params : the parameter locations for which the gradients are being
% computed. 
% RETURN g : gradients of the prior with respect to the parameters.
%
% SEEALSO : priorCreate
%
% COPYRIGHT : Neil D. Lawrence, 2003

% PRIOR
  
fhandle = str2func([prior.type 'PriorGradient']);
g = fhandle(prior, params);
