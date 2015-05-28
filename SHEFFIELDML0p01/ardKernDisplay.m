function ardKernDisplay(kern, spacing)

% ARDKERNDISPLAY Display parameters of the ARD kernel.
%
%	Description:
%
%	ARDKERNDISPLAY(KERN) displays the parameters of the pre-built RBF
%	and linear ARD kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	ARDKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	ARDKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2006 Neil D. Lawrence


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('RBF Variance: %2.4f\n', kern.rbfVariance)
fprintf(spacing);
fprintf('RBF inverse width: %2.4f\n', kern.inverseWidth)
fprintf(spacing);
fprintf('Linear Variance: %2.4f\n', kern.linearVariance)
fprintf(spacing);
fprintf('Noise Variance: %2.4f\n', kern.whiteVariance)
fprintf(spacing);
fprintf('Bias Variance: %2.4f\n', kern.biasVariance)
for i = 1:kern.inputDimension
  fprintf(spacing);
  fprintf('Input %d scale: %2.4f\n', i, kern.inputScales(i))
end
