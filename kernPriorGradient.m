function g = kernPriorGradient(kern)

% KERNPRIORGRADIENT Compute gradient terms associated with kernel priors.

% KERN

g = zeros(1, kern.nParams);
switch kern.type
 case {'cmpnd', 'multi', 'tensor'}
  startVal = 1;
  endVal = 0;
  for i = 1:length(kern.comp)
   endVal = endVal + kern.comp{i}.nParams;
   g(1, startVal:endVal) = kernPriorGradient(kern.comp{i});
   startVal = endVal + 1;
  end
  g = g*kern.paramGroups;
 otherwise
  if isfield(kern, 'priors')
    fhandle = str2func([kern.type 'KernExtractParam']);
    params = fhandle(kern);
    if iscell(kern.priors),
      for i = 1:length(kern.priors)
	index = kern.priors{i}.index;
	g(index) = g(index) + priorGradient(kern.priors{i}, params(index));
      end
    else
      for i = 1:length(kern.priors)
	index = kern.priors(i).index;
	g(index) = g(index) + priorGradient(kern.priors(i), params(index));
      end
    end

    factors = kernFactors(kern, 'gradfact');
    g(factors.index) = g(factors.index).*factors.val;
  end
end



  
