function [g, nu] = ncnmNoiseNuG(noise, mu, varSigma, y)

% NCNMNOISENUG Update nu and g parameters associated with null category noise model.
% FORMAT
% DESC computes the values of nu and g for use in IVM style updated
% equations. These are also used in EP style update equations. This
% command just calls ncnmNoiseGradVals, but then ensures that there
% are no nu values less than 0.
% ARG noise : the noise model structure.
% ARG mu : the input means to the noise model.
% ARG sigma : the input variances to the noise model.
% ARG y : the targets for the noise model.
% RETURN nu : the updated value for nu.
% RETURN g : the updated value for g.
%
% SEEALSO : ncnmNoiseGradVals, noiseUpdateNuG, ncnmNoiseParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE

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
