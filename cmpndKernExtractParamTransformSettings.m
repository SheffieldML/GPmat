function [paramtransformsettings, names] = cmpndKernExtractParamTransformSettings(kern)

% CMPNDKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings 
% from the CMPND kernel structure.
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
%
% COPYRIGHT : Antti Honkela, 2012

% KERN



paramtransformsettings=cell(size(kern.paramGroups,1),1);
if nargout > 1
  namesTemp = cell(1, size(kern.paramGroups,1));
end


startVal = 1;
endVal = 0;
storedTypes = cell(0);
for i = 1:length(kern.comp)
  [tempparamtransformsettings, tempnames] = kernExtractParamTransformSettings(kern.comp{i});
  endVal = endVal + length(tempparamtransformsettings);
  paramtransformsettings(startVal:endVal) = tempparamtransformsettings(:);

  if nargout > 1
    instNum = sum(strcmp(kern.comp{i}.type, storedTypes)) + 1;
    for j = startVal:endVal
      namesTemp{1, j} = [kern.comp{i}.type ' ' num2str(instNum) ' ' ...
                         tempnames{j - startVal + 1}];
    end;
    storedTypes{end+1} = kern.comp{i}.type;
  end;

  startVal = endVal + 1;
end


% If any parameters are 'tied together' deal with them.
paramtransformationsettingsFinal={};
paramGroups = kern.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  paramtransformationsettingsFinal{i} = paramtransformsettings{ind(1)};
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
paramtransformsettings = paramtransformationsettingsFinal;

