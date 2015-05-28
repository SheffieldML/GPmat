function kern = cmpndKernSetIndex(kern, component, indices)

% CMPNDKERNSETINDEX Set the indices in the compound kernel.
%
%	Description:
%
%	KERN = CMPNDKERNSETINDEX(KERN, COMPONENT, INDICES) sets the indices
%	of the input matrix to be used in the computation of the covariance
%	function.
%	 Returns:
%	  KERN - the covariance function with the indices set.
%	 Arguments:
%	  KERN - kernel matrix for which indices are to be set.
%	  COMPONENT - component number in the compound covariance.
%	  INDICES - indices of the input to be used.
%	
%
%	See also
%	KERNSETINDEX


%	Copyright (c) 2004, 2005, 2011 Neil D. Lawrence


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
