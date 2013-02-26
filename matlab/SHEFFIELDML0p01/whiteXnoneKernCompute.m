function K = whiteXnoneKernCompute(whiteKern, noneKern, x1, x2)

% WHITEXNONEKERNCOMPUTE Compute a cross kernel between WHITE and NONE kernels.
%
%	Description:
%
%	K = WHITEXNONEKERNCOMPUTE(WHITEKERN, NONEKERN, X) computes cross
%	kernel terms between white noise kernel and a dummy kernel for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white noise
%	   kernel.
%	  NONEKERN - the kernel structure associated with the dummy kernel.
%	  X - inputs for which kernel is to be computed.
%
%	K = WHITEXNONEKERNCOMPUTE(WHITEKERN, NONEKERN, X1, X2) computes
%	cross kernel terms between two WHITE kernels for the multiple output
%	kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white noise
%	   kernel.
%	  NONEKERN - the kernel structure associated with the dummy kernel.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, WHITEKERNPARAMINIT


%	Copyright (c) 2008 Neil D. Lawrence


if nargin < 4
  x2 = x1;
end

K = zeros(size(x1, 1), size(x2, 1));