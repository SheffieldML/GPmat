function K = whiteXwhiteKernCompute(whiteKern1, whiteKern2, x1, x2)

% WHITEXWHITEKERNCOMPUTE Compute a cross kernel between two WHITE kernels.
%
%	Description:
%
%	K = WHITEXWHITEKERNCOMPUTE(WHITEKERN1, WHITEKERN2, X) computes cross
%	kernel terms between two white noise kernels for the multiple output
%	kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  WHITEKERN1 - the kernel structure associated with the first white
%	   noise kernel.
%	  WHITEKERN2 - the kernel structure associated with the second white
%	   noise kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = WHITEXWHITEKERNCOMPUTE(WHITEKERN1, WHITEKERN2, X1, X2) computes
%	cross kernel terms between two WHITE kernels for the multiple output
%	kernel.
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
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, WHITEKERNPARAMINIT


%	Copyright (c) 2006 Neil D. Lawrence


if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));