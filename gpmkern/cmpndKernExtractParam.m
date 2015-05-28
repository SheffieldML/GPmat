function [params, names] = cmpndKernExtractParam(kern)

% CMPNDKERNEXTRACTPARAM Extract parameters from the CMPND kernel structure.
% FORMAT
% DESC Extract parameters from the compound kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of 'transforms' is assumed to be empty
% here. Any transformations of parameters should be done in
% component kernels.
%
% SEEALSO cmpndKernParamInit, cmpndKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


params = zeros(1, kern.nParams);
if nargout > 1
  namesTemp = cell(1, kern.nParams);
end
startVal = 1;
endVal = 0;
storedTypes = cell(0);
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if nargout > 1
    [params(1, startVal:endVal), namesTemp(startVal:endVal)] = kernExtractParam(kern.comp{i});
    instNum = sum(strcmp(kern.comp{i}.type, storedTypes)) + 1;
    for j = startVal:endVal
      namesTemp{1, j} = [kern.comp{i}.type ' ' num2str(instNum) ' ' ...
                         namesTemp{1, j}];
    end
    storedTypes{end+1} = kern.comp{i}.type;
  else
    params(1, startVal:endVal) = kernExtractParam(kern.comp{i});
  end
    startVal = endVal + 1;
end

% If any parameters are 'tied together' deal with them.
paramGroups = kern.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  if nargout > 1
    names{i} = namesTemp{ind(1)};
    if length(ind) > 1
      for j = 2:length(ind)
        names{i} = [names{i} ', ' namesTemp{ind(j)}];
      end
    end
  end
  paramGroups(ind(2:end), i) = 0;
end
params = params*paramGroups;
