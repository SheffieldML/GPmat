function [g, nu] = ncnmNoiseNuG(noise, mu, varSigma, y)

% NCNMNOISENUG Update nu and g parameters associated with null category noise model.
%
%	Description:
%
%	[NU, G] = NCNMNOISENUG(NOISE, MU, SIGMA, Y) computes the values of
%	nu and g for use in IVM style updated equations. These are also used
%	in EP style update equations. This command just calls
%	ncnmNoiseGradVals, but then ensures that there are no nu values less
%	than 0.
%	 Returns:
%	  NU - the updated value for nu.
%	  G - the updated value for g.
%	 Arguments:
%	  NOISE - the noise model structure.
%	  MU - the input means to the noise model.
%	  SIGMA - the input variances to the noise model.
%	  Y - the targets for the noise model.
%	
%
%	See also
%	NCNMNOISEGRADVALS, NOISEUPDATENUG, NCNMNOISEPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


[g, dlnZ_dvs] = feval([noise.type 'NoiseGradVals'], ...
		      noise, ...
		      mu, varSigma, ...
		      y);

nu = g.*g - 2*dlnZ_dvs;

% Reset any negative nu values to eps.
for i = 1:size(mu, 2)
  index = find(nu(:, i)< eps);
  nu(index) = eps;
end
