function [u, v] = kernPca(kern, X, numEig)

% KERNPCA performs KPCA.
% FORMAT
% DESC Performs kernel PCA on a given data set Y.
% ARG kern : the kernel type for performing kernel PCA.
% ARG Y : the data matrix with N rows and d columns.
% ARG dims : the number of dimensions of latent positions to
% return.
% RETURN U : the eigenvectors of the kernel matrix.
% RETURN lambda : the eigenvalues of the kernel matrix.
% 
% SEEALSO : kernCreate
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2006

% KERN

K = kernCompute(kern, X);
[u, v] = eigs(K, numEig);
