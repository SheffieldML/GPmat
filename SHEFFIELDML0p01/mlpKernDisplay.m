function mlpKernDisplay(kern, spacing)

% MLPKERNDISPLAY Display parameters of the MLP kernel.
%
%	Description:
%
%	MLPKERNDISPLAY(KERN) displays the parameters of the multi-layer
%	perceptron kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	MLPKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	MLPKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP kernel Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('MLP weight variance: %2.4f\n', kern.weightVariance)
fprintf(spacing);
fprintf('MLP bias variance: %2.4f\n', kern.biasVariance)