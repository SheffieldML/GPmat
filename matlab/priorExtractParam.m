function [params, names] = priorExtractParam(prior)

% PRIOREXTRACTPARAM Extract the prior model's parameters.

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