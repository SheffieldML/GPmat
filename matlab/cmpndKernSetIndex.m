function kern = cmpndKernSetIndex(kern, component, indices)

% CMPNDKERNSETINDEX Set the indices in the compound kernel.

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
kern.comp{component} = kernParamInit(kern.comp{component});
