function invcmpndKernDisplay(kern, varargin)

% INVCMPNDKERNDISPLAY Display parameters of the INVCMPND kernel.
% FORMAT
% DESC displays the parameters of the inv. precisions compound
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : invcmpndKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, Andreas C. Damianou, 2012

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
fprintf('Inverse Precision Matrix Compound kernel:\n')
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end
