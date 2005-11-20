function [g, nu] = noiseUpdateNuG(noise, mu, varSigma, y);

% NOISEUPDATENUG Update nu and g for a given noise model.

% NOISE

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

