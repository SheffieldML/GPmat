
% PCAEMBED Embed data set with PCA.
%
%	Description:
%
%	[X, U, V] = PCAEMBED(Y, DIMS)  [X, U, V] = PCAEMBED(DATA, N) returns
%	latent positions for a given data set via PCA as well as principal
%	components and eigenvalues.
%	 Returns:
%	  X - the latent positions.
%	  U - the largest N eigenvalues
%	  V - the principal components
%	 Arguments:
%	  Y - the data set which you want the latent positions for.
%	  DIMS - the dimensionality of the latent space.
%	
%
%	See also
%	PCA


%	Copyright (c) 2012 Andreas Damianou


function [X, U, V] = pcaEmbed(Y, dims)
    [U,V] = pca(Y,dims);
    X = Y*V;
