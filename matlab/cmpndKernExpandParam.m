function kern = cmpndKernExpandParam(kern, params)

% CMPNDKERNEXPANDPARAM Create kernel structure from parameters.

% KERN

params = params*kern.paramGroups';
startVal = 1;
endVal = 0;
kern.whiteVariance = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  kern.comp{i} = kernExpandParam(kern.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
  if strcmp(kern.comp{i}.type, 'white')
    kern.whiteVariance = kern.whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      kern.whiteVariance = kern.whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end
