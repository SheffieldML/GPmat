function ivmDisplay(model, spacing)

% IVMDISPLAY Display parameters of an IVM model.
% FORMAT
% DESC displays the parameters of an IVM model.
% ARG model : the model to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG model : the IVM model to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : gaussianNoiseParamInit, modelDisplay, noiseDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2007

% IVM

if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('IVM Model:\n')
fprintf(spacing);
fprintf(' Noise Model:\n')
noiseDisplay(model.noise, length(spacing)+2);
fprintf(spacing);
fprintf(' Kernel:\n');
kernDisplay(model.kern, length(spacing)+2);
