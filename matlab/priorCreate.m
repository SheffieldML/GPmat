function prior = priorCreate(type)

% PRIORCREATE Create a prior structure given a type.
% FORMAT
% DESC creates a prior structure given a type.
% ARG type : Type of prior to be created,  some standard types are
% 'gamma', 'gaussian', 'laplace' and 'invgamma'.
% RETURN prior : The prior structure.
%
% SEEALSO : priorParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2003

% SHEFFIELDML

prior.type = type;
prior = priorParamInit(prior);
