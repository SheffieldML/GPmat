function [params, names] = priorExtractParam(prior)

% PRIOREXTRACTPARAM Extract the prior model's parameters.
% FORMAT
% DESC Extract parameters from the prior into a vector of
% parameters for optimisation.
% ARG prior : the prior structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the prior. If
% the field 'transforms' is not empty in the prior, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO : priorExpandParam, scg, conjgrad
% 
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005

% PRIOR

fhandle = str2func([prior.type 'PriorExtractParam']);
if nargout < 2
  params = fhandle(prior);
else
  [params, names] = fhandle(prior);
end


% Check if parameters are being optimised in a transformed space.
if isfield(prior, 'transforms')
  for i = 1:length(prior.transforms)
    index = prior.transforms(i).index;
    fhandle = str2func([prior.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'xtoa');
  end
end
