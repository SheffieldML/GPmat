function translateKernDisplay(kern, varargin)

% TRANSLATEKERNDISPLAY Display parameters of the TRANSLATE kernel.
% FORMAT
% DESC displays the parameters of the input space translation
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : translateKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

if nargin > 1
  spacing = repmat(32, 1, varargin{1});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Translate kernel:\n')
for i = 1:length(kern.centre)
  fprintf(spacing);
  fprintf(' Centre %d: %2.4f\n', i, kern.centre(i));
end
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end
