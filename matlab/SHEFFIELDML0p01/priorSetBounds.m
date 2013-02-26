function prior = priorSetBounds(prior, bounds)

% PRIORSETBOUNDS Set the bounded prior model's bounds from bounds vector.
%
%	Description:
%
%	PRIOR = PRIORSETBOUNDS(PRIOR, BOUNDS) sets the bounds of a bounded
%	prior.
%	 Returns:
%	  PRIOR - prior structure with the given bounds in the relevant
%	   locations.
%	 Arguments:
%	  PRIOR - the prior structure in which the bounds are to be placed.
%	  BOUNDS - vector of bounds which are to be placed in the prior
%	   structure.
%	
%
%	See also
%	PRIOREXPANDPARAM


%	Copyright (c) 2011 Antti Honkela


if iscell(bounds) && numel(bounds)==1,
  bounds = bounds{1};
end

if isfield(prior, 'isBounded') && prior.isBounded,
  fhandle = str2func([prior.type 'PriorSetBounds']);
  prior = fhandle(prior, bounds);
end
