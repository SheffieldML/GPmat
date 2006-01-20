function cmpndKernDisplay(kern, varargin)

% CMPNDKERNDISPLAY Display the parameters of the compound kernel.

% KERN

for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i}, varargin{:});
end