% PCAEMBED Embed data set with PCA.
% FORMAT
% DESC 
% [X, U, V] = PCAEMBED(DATA, N) returns latent positions for a given data set via PCA
% as well as principal components and eigenvalues. 
% ARG Y : the data set which you want the latent positions for.
% ARG dims : the dimensionality of the latent space.
% RETURN X : the latent positions.
% RETURN U : the largest N eigenvalues
% RETURN V : the principal components
% 
% COPYRIGHT : Andreas Damianou, 2012
%
% SEEALSO : pca

% MLTOOLS

function [X, U, V] = pcaEmbed(Y, dims)
    [U,V] = pca(Y,dims);
    X = Y*V;
