function expKernDisplay(kern, spacing)

% EXPKERNDISPLAY Display parameters of the EXP kernel.
% FORMAT
% DESC displays the parameters of the exponentiated
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : expKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

if nargin > 1
  spacing = repmat(32, 1, varargin{:});
  varargin{1} = varargin{1}+2;
else
  spacing = [];
  varargin{1} = 2;
end
spacing = char(spacing);
fprintf(spacing);
fprintf('Exponentiated kernel:\n')
fprintf(spacing);
fprintf('EXP variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('Log process kernels:\n')
kernDisplay(kern.argument, varargin{:});
