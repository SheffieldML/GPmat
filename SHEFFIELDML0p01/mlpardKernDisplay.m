function mlpardKernDisplay(kern, spacing)

% MLPARDKERNDISPLAY Display parameters of the MLPARD kernel.
%
%	Description:
%
%	MLPARDKERNDISPLAY(KERN) displays the parameters of the automatic
%	relevance determination multi-layer perceptron kernel and the kernel
%	type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	MLPARDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	MLPARDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence



if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD kernel Variance: %2.4f\n', kern.variance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD weight variance: %2.4f\n', kern.weightVariance)
spacing = char(spacing);
fprintf(spacing);
fprintf('MLP ARD bias variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  spacing = char(spacing);
  fprintf(spacing);
  fprintf('MLP ARD Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
