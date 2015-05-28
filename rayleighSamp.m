function y = rayleighSamp(sigma, numSamps)

% RAYLEIGHSAMP Sample from a Rayleigh with a given sigma.
% FORMAT 
% DESC samples a given number of samples from a Rayleigh with a
% given sigma matrix.
% ARG sigma : the sigma of the Rayleigh to sample from.
% ARG numSamps : the number of samples to take from Rayleigh.
% RETURN y : the samples from the Rayleigh
%
% SEEALSO : rand, gaussSamp
%
% COPYRIGHT : Neil D. Lawrence, 2012

% NDLUTIL

y = rand(numSamps, 1);
y = sigma.*sqrt(-2*log(y));
