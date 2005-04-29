function y = noiseOut(noise, mu, varsigma);

% NOISEOUT Give the output of the noise model given the mean and variance.

% NOISE

fhandle = str2func([noise.type 'NoiseOut']);
if str2num(version('-release'))>13
  y = fhandle(noise, mu, varsigma);
else
  y = feval(fhandle, noise, mu, varsigma);
end