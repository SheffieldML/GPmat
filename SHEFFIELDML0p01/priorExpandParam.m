function prior = priorExpandParam(prior, params)

% PRIOREXPANDPARAM Expand the prior model's parameters from params vector.
%
%	Description:
%
%	PRIOR = PRIOREXPANDPARAM(PRIOR, PARAM) returns a prior structure
%	filled with the parameters in the given vector. This is used as a
%	helper function to enable parameters to be optimised in, for
%	example, the NETLAB optimisation functions.
%	 Returns:
%	  PRIOR - prior structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  PRIOR - the prior structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the prior
%	   structure.
%	As well as extracting the parameters, some transformation of
%	parameters is also undertaken in this file. If the field
%	transforms is not empty, it dictactes how the prior parameters
%	are to be transformed (for example by a exponential to keep them
%	positive).
%	
%	
%
%	See also
%	PRIOREXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2003, 2004 Neil D. Lawrence



if isfield(prior, 'transforms')
  for i = 1:length(prior.transforms)
    index = prior.transforms(i).index;
    fhandle = str2func([prior.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end

fhandle = str2func([prior.type 'PriorExpandParam']);
prior = fhandle(prior, params);
