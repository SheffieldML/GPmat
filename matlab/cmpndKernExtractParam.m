function params = cmpndKernExtractParam(kern)


% CMPNDKERNEXTRACTPARAM Extract parameters from the CMPND kernel structure.
% FORMAT
% DESC Extract parameters from the compound kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO cmpndKernParamInit, cmpndKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% KERN


params = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  params(1, startVal:endVal) = kernExtractParam(kern.comp{i});
  startVal = endVal + 1;
end
paramGroups = kern.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  paramGroups(ind(2:end), i) = 0;
end
params = params*paramGroups;
