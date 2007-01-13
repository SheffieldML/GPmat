function orderedNoiseDisplay(noise, spacing)


% ORDEREDNOISEDISPLAY Display parameters of the ORDERED noise.
% FORMAT
% DESC displays the parameters of the ordered categorical
% noise and the noise type to the console.
% ARG noise : the noise to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG noise : the noise to display.
% ARG spacing : how many spaces to indent the display of the noise by.
%
% SEEALSO : orderedNoiseParamInit, modelDisplay, noiseDisplay
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
  fprintf('Ordered bias on process %d: %2.4f\n', i, noise.bias(i))
end
for i = 1:noise.C-2
  fprintf('Ordered noise model width %d: %2.4f\n', i, noise.widths(i))
end
fprintf(spacing);
fprintf(spacing);
fprintf('Ordered Sigma2: %2.4f\n', noise.variance);