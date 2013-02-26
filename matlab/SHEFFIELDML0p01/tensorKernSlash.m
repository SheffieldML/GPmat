function kern = tensorKernSlash(kern, index)

% TENSORKERNSLASH Tensor kernel created by removing ith component.
%
%	Description:
%	kern = tensorKernSlash(kern, index)
%

kern.nParams = kern.nParams - kern.comp{index}.nParams;
kern.comp = kern.comp([1:index-1 index+1:end]);
