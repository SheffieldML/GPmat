function k = rbfhKernDiagCompute(kern, x)

% RBFHKERNDIAGCOMPUTE Compute diagonal of RBFH kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the radial basis 
% function heat kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : rbfhKernParamInit, kernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Since the variance is on the diagonal is just ones

if size(x, 2) ~= 2
    error('Input can only have two columns');
end

% Split the domain into time domain and spatial domain and account for
% missing values. If there are no missing values the computation of the
% kernel is a pointwise prodruct, otherwise it is a kronecker product.
t = x(x(:,1)~=Inf,1);
s = x(x(:,2)~=Inf,2);

if (length(t) == length(s))
    ut = unique(t);
    us = unique(s);
    if (length(ut)*length(us) == length(t))
        k = ones(length(ut)*length(us),1);        
    else
        k = ones(length(t),1);        
    end    
else
    k = ones(length(t)*length(s),1);
end
