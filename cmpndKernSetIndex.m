function kern = cmpndKernSetIndex(kern, component, indices)

% CMPNDKERNSETINDEX Set the indices in the compound kernel.
% FORMAT
% DESC sets the indices of the input matrix to be used in the computation
% of the covariance function.
% ARG kern : kernel matrix for which indices are to be set.
% ARG component : component number in the compound covariance.
% ARG indices : indices of the input to be used.
% RETURN kern : the covariance function with the indices set.
% 
% SEEALSO : kernSetIndex
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2011

% KERN

if size(indices, 1) ~= 1
  error('Indices should be a row vector.');
end
if max(indices) > kern.inputDimension
  error('Indices are larger than kernel input dimension');
elseif min(indices) < 0
  error('Indices should be greater than zero.');
end
  
kern.comp{component}.inputDimension = length(indices);
kern.comp{component}.index = indices;
param = kernExtractParam(kern.comp{component});
kern.comp{component} = kernParamInit(kern.comp{component});
kern.comp{component} = kernExpandParam(kern.comp{component}, param);

kern.nParams = 0;
for i = 1:length(kern.comp)
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
end


kern.paramGroups = speye(kern.nParams);
