function [g, nu] = noiseUpdateNuG(noise, mu, varSigma, y);

% NOISEUPDATENUG Update nu and g for a given noise model.
%
%	Description:
%
%	[G, NU] = NOISEUPDATENUG(NOISE, MU, VARSIGMA, Y, Y) computes the
%	values nu and g for the given noise given the mean and variance
%	inputs as well as the output of the noise model.
%	 Returns:
%	  G - the vector g, which is the gradient of log Z with respect to
%	   the input mean.
%	  NU - the vector nu, see equation 10 of "Extensions of the
%	   Informative Vector Machine".
%	 Arguments:
%	  NOISE - the noise structure for which the nu and g are computed.
%	  MU - input mean to the noise model.
%	  VARSIGMA - input variance to the noise model.
%	  Y - target output for the noise model.
%	  Y - target output for the noise model.
%	
%
%	See also
%	NOISEPARAMINIT, NOISEUPDATESITES, NOISECREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence


if noise.updateNuG
  % The noise model has it's own code for site updates.
  fhandle = str2func([noise.type 'NoiseNuG']);
  [g, nu] = fhandle(noise, mu, varSigma, y);
else
  % Use the standard (general) code.
  fhandle = str2func([noise.type 'NoiseGradVals']);
  [g, dlnZ_dvs] = fhandle(noise, mu, varSigma, y);
  nu = g.*g - 2*dlnZ_dvs;
  nu(find(abs(nu) < eps)) = eps;
end

