function X = smallrandEmbed(Y, dims)

% SMALLRANDEMBED Embed data set with small random values.
% FORMAT
% DESC returns (initial) latent positions for a given data set.
% ARG Y : the data set which you want the latent positions for.
% ARG dims : the dimensionality of the latent space.
% RETURN X : the latent positions.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% SEEALSO : lleEmbed, isomapEmbed, ppcaEmbed

% MLTOOLS

X = randn(size(Y, 1), dims)*0.0001;
