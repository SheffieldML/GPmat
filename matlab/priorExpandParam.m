function prior = priorExpandParam(prior, params)

% PRIOREXPANDPARAM Expand the prior model's parameters from params vector.
% IVM

if isfield(prior, 'transforms')
  for i = 1:length(prior.transforms)
    index = prior.transforms(i).index;
    params(index) = feval([prior.transforms(i).type 'Transform'], ...
              params(index), 'atox');
  end
end

prior = feval([prior.type 'PriorExpandParam'], prior, params);
