function [m, beta] = noiseUpdateSites(noise, g, nu, mu, varSigma, y);

% NOISEUPDATESITES Update site parameters for a given noise model.
% FORMAT
% DESC updates the site parameters for the given
% noise model. 
% ARG noise : the noise structure for which the site parameters are to
% be updated. 
% ARG g : values of g as retuned by noiseUpdateNuG.
% ARG nu : values of nu as retuned by noiseUpdateNuG.
% ARG mu : the mean value of the Gaussian input to the noise structure.
% ARG varSigma : the variance of the Gaussian input to the noise structure.
% ARG y : the target value.
% RETURN m : the site mean parameters.
% RETURN beta : the site precision parameters.
%
% SEEALSO : noiseParamInit, noiseUpdateNuG
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

if noise.updateSites
  % The noise model has its own code for site updates.
  fhandle = str2func([noise.type 'NoiseSites']);
  [m, beta] = fhandle(noise, g, nu,  mu, varSigma, y);
else
  % Use the standard code.
  beta = nu./(1-nu.*varSigma);
  m = mu + g./nu;
end
