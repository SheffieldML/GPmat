function [m, beta] = noiseUpdateSites(noise, g, nu, mu, varSigma, y);

% NOISEUPDATESITES Update site parameters for a given noise model.
%
%	Description:
%
%	[M, BETA] = NOISEUPDATESITES(NOISE, G, NU, MU, VARSIGMA, Y) updates
%	the site parameters for the given noise model.
%	 Returns:
%	  M - the site mean parameters.
%	  BETA - the site precision parameters.
%	 Arguments:
%	  NOISE - the noise structure for which the site parameters are to
%	   be updated.
%	  G - values of g as retuned by noiseUpdateNuG.
%	  NU - values of nu as retuned by noiseUpdateNuG.
%	  MU - the mean value of the Gaussian input to the noise structure.
%	  VARSIGMA - the variance of the Gaussian input to the noise
%	   structure.
%	  Y - the target value.
%	
%
%	See also
%	NOISEPARAMINIT, NOISEUPDATENUG


%	Copyright (c) 2004, 2005 Neil D. Lawrence


if noise.updateSites
  % The noise model has its own code for site updates.
  fhandle = str2func([noise.type 'NoiseSites']);
  [m, beta] = fhandle(noise, g, nu,  mu, varSigma, y);
else
  % Use the standard code.
  beta = nu./(1-nu.*varSigma);
  m = mu + g./nu;
end
