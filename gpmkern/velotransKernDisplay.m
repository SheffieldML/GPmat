function velotransKernDisplay(kern, varargin)

% VELOTRANSKERNDISPLAY Display parameters of the VELOTRANS kernel.
% FORMAT
% DESC displays the parameters of the velocity translate
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : velotransKernParamInit, modelDisplay, kernDisplay, translateKernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2011

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
fprintf('Velocity Translate kernel:\n')
for i = 1:length(kern.velocity)
  fprintf(spacing);
  fprintf(' Velocity %d: %2.4f\n', i, kern.velocity(i));
end
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end

  
