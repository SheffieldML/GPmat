function sparseKernDisplay(kern, varargin)

% SPARSEKERNDISPLAY Display parameters of the SPARSE kernel.
% FORMAT
% DESC displays the parameters of the sparse output block
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : sparseKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

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
fprintf('Sparse kernel:\n')
for i = 1:length(kern.comp)
  fprintf(spacing);
  fprintf('Block %d\n', i)
  kernDisplay(kern.comp{i}, varargin{:});
end
