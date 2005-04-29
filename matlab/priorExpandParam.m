function prior = priorExpandParam(prior, params)

% PRIOREXPANDPARAM Expand the prior model's parameters from params vector.

% PRIOR


if isfield(prior, 'transforms')
  for i = 1:length(prior.transforms)
    index = prior.transforms(i).index;
    fhandle = str2func([prior.transforms(i).type 'Transform']);
    params(index) = fhandle(params(index), 'atox');
  end
end

fhandle = str2func([prior.type 'PriorExpandParam']);
prior = fhandle(prior, params);
