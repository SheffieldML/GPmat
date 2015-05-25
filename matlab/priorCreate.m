function prior = priorCreate(type, options)

% PRIORCREATE Create a prior structure given a type.
% FORMAT
% DESC creates a prior structure given a type.
% ARG type : Type of prior to be created,  some standard types are
% 'gamma', 'gaussian', 'laplace' and 'invgamma'.
% ARG options : Options for the prior (type specific).
% RETURN prior : The prior structure.
%
% SEEALSO : priorParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2003
%
% COPYRIGHT : Antti Honkela, 2013

% PRIOR

prior.type = type;
if nargin < 2,
  prior = priorParamInit(prior);
else
  prior = priorParamInit(prior, options);
end
