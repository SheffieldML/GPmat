function [k, n2] = fileKernCompute(kern, varargin)

% FILEKERNCOMPUTE Compute the kernel given the parameters and indices.

% KERN

k = kern.variance*fileKernRead(kern, varargin{:});
