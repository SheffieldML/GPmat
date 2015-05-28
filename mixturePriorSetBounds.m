function prior = mixturePriorSetBounds(prior, bounds)

% MIXTUREPRIORSETBOUNDS Set mixture prior bounds.

% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

if prior.isBounded,
  for k=1:length(prior.comp),
    prior.comp{k} = priorSetBounds(prior.comp{k}, bounds);
  end
else
  warning('Trying to set bounds on unbounded prior ignored.')
end
