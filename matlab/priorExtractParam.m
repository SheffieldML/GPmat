function [params, names] = priorExtractParam(prior)

% PRIOREXTRACTPARAM Extract the prior model's parameters.

% PRIOR

% PRIOR


if nargout < 2
  params = feval([prior.type 'PriorExtractParam'], prior);
else
  [params, names] = feval([prior.type 'PriorExtractParam'], prior);
end


% Check if parameters are being optimised in a transformed space.
if isfield(prior, 'transforms')
  for i = 1:length(prior.transforms)
    index = prior.transforms(i).index;
    params(index) = feval([prior.transforms(i).type 'Transform'], ...
              params(index), 'xtoa');
  end
end