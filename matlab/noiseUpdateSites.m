function [m, beta] = noiseUpdateSites(noise, g, nu, mu, varSigma, y);

% NOISEUPDATESITES Update site parameters for a given noise model.

% NOISE

if noise.updateSites
  % The noise model has it's own code for site updates.
  fhandle = str2func([noise.type 'NoiseSites']);
  [m, beta] = fhandle(noise, g, nu,  mu, varSigma, y);
else
  % Use the standard code.
  beta = nu./(1-nu.*varSigma);
  m = mu + g./nu;
end
