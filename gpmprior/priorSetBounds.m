function prior = priorSetBounds(prior, bounds)

% PRIORSETBOUNDS Set the bounded prior model's bounds from bounds vector.
% FORMAT
% DESC sets the bounds of a bounded prior.
% ARG prior : the prior structure in which the bounds are to be
% placed.
% ARG bounds : vector of bounds which are to be placed in the
% prior structure.
% RETURN prior : prior structure with the given bounds in the
% relevant locations.
%
% SEEALSO : priorExpandParam
%
% COPYRIGHT : Antti Honkela, 2011

% PRIOR

if iscell(bounds) && numel(bounds)==1,
  bounds = bounds{1};
end

if isfield(prior, 'isBounded') && prior.isBounded,
  fhandle = str2func([prior.type 'PriorSetBounds']);
  prior = fhandle(prior, bounds);
end
