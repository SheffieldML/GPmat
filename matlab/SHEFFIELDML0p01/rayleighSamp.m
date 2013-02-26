function y = rayleighSamp(sigma, numSamps)

% RAYLEIGHSAMP Sample from a Rayleigh with a given sigma.
%
%	Description:
%
%	Y = RAYLEIGHSAMP(SIGMA, NUMSAMPS) samples a given number of samples
%	from a Rayleigh with a given sigma matrix.
%	 Returns:
%	  Y - the samples from the Rayleigh
%	 Arguments:
%	  SIGMA - the sigma of the Rayleigh to sample from.
%	  NUMSAMPS - the number of samples to take from Rayleigh.
%	
%
%	See also
%	RAND, GAUSSSAMP


%	Copyright (c) 2012 Neil D. Lawrence


y = rand(numSamps, 1);
y = sigma.*sqrt(-2*log(y));
