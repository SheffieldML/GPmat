function prior = gaussianPriorExpandParam(prior, params)

% GAUSSIANPRIOREXPANDPARAM Expand Gaussian prior structure from param vector.
% FORMAT
% DESC places the given parameters in a Gaussian prior structure.
% ARG prior : the structure to place the parameters in.
% ARG params : the parameters to place in the structure.
% RETURN prior : the structure with the parameters in place.
%
% SEEALSO priorExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004

% PRIOR

prior.precision = params(1);
