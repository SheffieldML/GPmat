function L = kernPriorLogProb(kern)

% KERNPRIORLOGPROB Compute penalty terms associated with kernel priors.

% KERN

L = 0;
switch kern.type
 case 'cmpnd'
  for i = 1:length(kern.comp)
    L = L + kernPriorLogProb(kern.comp{i});
  end
 otherwise
  if isfield(kern, 'priors')
    params = feval([kern.type 'KernExtractParam'], kern);
    for i = 1:length(kern.priors)
      index = kern.priors(i).index;
      L = L + priorLogProb(kern.priors(i), params(index));
    end
  end
end