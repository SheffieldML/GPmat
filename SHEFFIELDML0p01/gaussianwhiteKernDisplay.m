function gaussianwhiteKernDisplay(kern, spacing)

% GAUSSIANWHITEKERNDISPLAY Display parameters of the GAUSSIAN white kernel.
%
%	Description:
%
%	GAUSSIANWHITEKERNDISPLAY(KERN) displays the parameters of the
%	gaussian white kernel and the kernel type to the console.
%	 Arguments:
%	  KERN - the kernel to display.
%
%	GAUSSIANWHITEKERNDISPLAY(KERN, SPACING)
%	 Arguments:
%	  KERN - the kernel to display.
%	  SPACING - how many spaces to indent the display of the kernel by.
%	
%	
%
%	See also
%	GAUSSIANWHITEKERNPARAMINIT, MODELDISPLAY, KERNDISPLAY


%	Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009.


if nargin > 1
  spacing = repmat(32, 1, spacing);
else
  spacing = [];
end
spacing = char(spacing);
fprintf(spacing);
if kern.isArd
    limit = kern.inputDimension;
else
    limit = 1;
end
for j =1:kern.nIndFunct
    for k =1:limit,
        fprintf('VIK %d inverse width %d : %2.4f (length scale %2.4f)\n', ...
            j, k, kern.precisionT(k,j), 1/sqrt(kern.precisionT(k,j)));
        fprintf(spacing);
    end
end
fprintf('White noise variance: %2.4f\n', kern.sigma2Noise)
