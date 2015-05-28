function probitNoiseDisplay(noise, spacing)


% PROBITNOISEDISPLAY Display parameters of the PROBIT noise.
% FORMAT
% DESC displays the parameters of the probit based classification
% noise and the noise type to the console.
% ARG noise : the noise to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG noise : the noise to display.
% ARG spacing : how many spaces to indent the display of the noise by.
%
% SEEALSO : probitNoiseParamInit, modelDisplay, noiseDisplay
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
  fprintf('Probit bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf(spacing);
fprintf('Probit Sigma2: %2.4f\n', noise.sigma2);


