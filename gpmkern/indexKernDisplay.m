function indexKernDisplay(kern, spacing)

% INDEXKERNDISPLAY Display parameters of the INDEX kernel.
% FORMAT
% DESC displays the parameters of the index based covariance function
% kernel and the kernel type to the console.
% ARG kern : the kernel to display.
%
% FORMAT does the same as above, but indents the display according
% to the amount specified.
% ARG kern : the kernel to display.
% ARG spacing : how many spaces to indent the display of the kernel by.
%
% SEEALSO : indexKernParamInit, modelDisplay, kernDisplay
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN


  if nargin > 1
    spacing = repmat(32, 1, spacing);
  else
    spacing = [];
  end
  spacing = char(spacing);
  fprintf(spacing)
  fprintf('Index kernel Variance: %2.4f\n', kern.variance)
end
