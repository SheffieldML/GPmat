function ncnmNoiseDisplay(noise, spacing)

% NCNMNOISEDISPLAY Display  parameters from null category noise model.
% FORMAT
% DESC displays the parameters of the null category noise model.
% ARG noise : the noise model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG noise : the noise model to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO ncnmNoiseParamInit, modelDisplay, noiseDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NOISE

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
for i = 1:noise.numProcess
  fprintf(spacing);
  fprintf('Null Cat bias on process %d: %2.4f\n', i, noise.bias(i))
end
fprintf(spacing);
fprintf('Null cat width: %2.4f\n', noise.width)
fprintf(spacing);
fprintf('Null cat Gamma -: %2.4f\n', noise.gamman)
fprintf(spacing);
fprintf('Null cat Gamma +: %2.4f\n', noise.gammap)
fprintf(spacing);
fprintf('Null cat Sigma2: %2.4f\n', noise.sigma2)
