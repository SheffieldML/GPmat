function disimKernDisplay(kern, spacing)

% DISIMKERNDISPLAY Display parameters of the DISIM kernel.
%
%	Description:
%
%	DISIMKERNDISPLAY(KERN) displays the parameters of the driven input
%	single input motif kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	DISIMKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	DISIMKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2006 Neil D. Lawrence
%	Copyright (c) 2007 Antti Honkela


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
fprintf('DISIM Variance: %2.4f\n', kern.di_variance)
fprintf(spacing);
fprintf('DISIM decay: %2.4f\n', kern.di_decay)
fprintf(spacing);
fprintf('SIM Variance: %2.4f\n', kern.variance)
fprintf(spacing);
fprintf('SIM inverse width: %2.4f (length scale %2.4f)\n', ...
        kern.inverseWidth, 1/sqrt(kern.inverseWidth));
fprintf(spacing);
fprintf('SIM decay: %2.4f\n', kern.decay)
fprintf(spacing);
fprintf('RBF variance: %2.4f\n', kern.rbf_variance)
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  fprintf(spacing);
  fprintf('SIM Initial Variance: %2.4f\n', kern.initialVariance)
end
