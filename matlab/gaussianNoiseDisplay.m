function gaussianNoiseDisplay(noise, spacing)

% GAUSSIANNOISEDISPLAY Display parameters of the GAUSSIAN noise.
% FORMAT
% DESC displays the parameters of the Gaussian
% noise and the noise type to the console.
% ARG noise : the noise to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG noise : the noise to display.
% ARG spacing : how many spaces to indent the display of the noise by.
%
% SEEALSO : gaussianNoiseParamInit, modelDisplay, noiseDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i = 1:noise.numProcess
  fprintf(spacing);
  fprintf('Gaussian bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf(spacing);
fprintf('Gaussian noise: %2.4f\n', noise.sigma2);