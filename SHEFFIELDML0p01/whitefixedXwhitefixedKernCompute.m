function K = whitefixedXwhitefixedKernCompute(whiteKern1, whiteKern2, x1, x2)

% WHITEFIXEDXWHITEFIXEDKERNCOMPUTE Compute a cross kernel between two WHITEFIXED kernels.
%
%	Description:
%
%	K = WHITEFIXEDXWHITEFIXEDKERNCOMPUTE(WHITEKERN1, WHITEKERN2, X)
%	computes cross kernel terms between two whitefixed noise kernels for
%	the multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  WHITEKERN1 - the kernel structure associated with the first white
%	   noise kernel.
%	  WHITEKERN2 - the kernel structure associated with the second white
%	   noise kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = WHITEFIXEDXWHITEFIXEDKERNCOMPUTE(WHITEKERN1, WHITEKERN2, X1, X2)
%	computes cross kernel terms between two WHITEFIXED kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  WHITEKERN1 - the kernel structure associated with the first white
%	   noise kernel.
%	  WHITEKERN2 - the kernel structure associated with the second white
%	   noise kernel.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, WHITEFIXEDKERNPARAMINIT


%	Copyright (c) 2008 Neil D. Lawrence


if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));