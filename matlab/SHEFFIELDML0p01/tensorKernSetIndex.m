function kern = tensorKernSetIndex(kern, component, indices)

% TENSORKERNSETINDEX Set the indices in the tensor kernel.
%
%	Description:
%	kern = tensorKernSetIndex(kern, component, indices)
%

kern = cmpndKernSetIndex(kern, component, indices);