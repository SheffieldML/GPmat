function tensorKernDisplay(kern, varargin)

% TENSORKERNDISPLAY Display the parameters of the compound kernel.

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
fprintf('Tensor kernel:\n')
for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end