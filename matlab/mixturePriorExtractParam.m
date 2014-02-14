function [params, names] = mixturePriorExtractParam(prior)

% MIXTUREPRIOREXTRACTPARAM Extract params from mixture prior structure.

% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

params = zeros(1, prior.nParams);
if nargout > 1
  names = cell(1, prior.nParams);
end

begindex = 1;
for k=1:length(prior.comp),
  myindex = begindex:(begindex+prior.comp{k}.nParams-1);
  if nargout > 1,
    [params(myindex), names(myindex)] = priorExtractParam(prior.comp{k});
  else
    params(myindex) = priorExtractParam(prior.comp{k});
  end
  begindex = begindex + prior.comp{k}.nParams;
end
