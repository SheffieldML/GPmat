function kern = invcmpndKernSetIndex(kern, component, indices)

% INVCMPNDKERNSETINDEX Set the indices in the inv. compound kernel.
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
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN


kern = cmpndKernSetIndex(kern, component, indices);
