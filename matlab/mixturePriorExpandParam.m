function prior = mixturePriorExpandParam(prior, params)

% MIXTUREPRIOREXPANDPARAM Expand mixture prior structure from params.

% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

begindex = 1;
for k=1:length(prior.comp),
  myindex = begindex:(begindex+prior.comp{k}.nParams-1);
  if prior.comp{k}.nParams > 0,
    prior.comp{k} = priorExpandParam(prior.comp{k}, params(myindex));
  end
  begindex = begindex + prior.comp{k}.nParams;
end
