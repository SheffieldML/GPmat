function sqexpKernDisplay(kern, spacing)

% SQEXPKERNDISPLAY Display parameters of the SQEXP kernel.
%
%	Description:
%
%	SQEXPKERNDISPLAY(KERN) displays the parameters of the pre-built
%	compound squared exponential kernel and the kernel type to the
%	console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	SQEXPKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%
%	See also
%	SQEXPKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004 Neil D. Lawrence



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
fprintf('White noise Variance: %2.4f\n', kern.whiteVariance)
fprintf(spacing);
fprintf('Bias Kernel Variance: %2.4f\n', kern.biasVariance)
