function indexardKernDisplay(kern, spacing)

% INDEXARDKERNDISPLAY Display parameters of the INDEXARD kernel.
% FORMAT
% DESC displays the parameters of the index ard based covariance function
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
  for i = 1:length(kern.indices)
    fprintf(spacing)
    fprintf('Index ARD Index %d scale: %2.4f\n', kern.indices(i), kern.indexScales(i))
  end
end
