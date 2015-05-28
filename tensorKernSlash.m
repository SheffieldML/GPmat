function kern = tensorKernSlash(kern, index)

% TENSORKERNSLASH Tensor kernel created by removing ith component.

% KERN

kern.nParams = kern.nParams - kern.comp{index}.nParams;
kern.comp = kern.comp([1:index-1 index+1:end]);
