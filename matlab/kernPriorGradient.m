function g = kernPriorGradient(kern)

% KERNPRIORGRADIENT Compute gradient terms associated with kernel priors.

% KERN

g = zeros(1, kern.nParams);
switch kern.type
 case 'cmpnd'
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
    for i = 1:length(kern.priors)
      index = kern.priors(i).index;
      g(index) = g(index) + priorGradient(kern.priors(i), params(index));
    end
    % Check if parameters are being optimised in a transformed space.
    if isfield(kern, 'transforms')
      for i = 1:length(kern.transforms)
        index = kern.transforms(i).index;
        fhandle = str2func([kern.transforms(i).type 'Transform']);
        g(index) = g(index).*fhandle(params(index), 'gradfact');
      end
    end
  end
end



  