function k = invcmpndKernCompute(kern, x, x2)

% INVCMPNDKERNCOMPUTE Compute the INVERSE-PRECISION-CMPND kernel given the parameters and X.
%
%	Description:
%
%	K = INVCMPNDKERNCOMPUTE(KERN, X, X2) computes the kernel parameters
%	for the kernel that is the inverse of the sum of the precisions of 1
%	or more kernels, given inputs associated with rows and columns. That
%	is, K(x,x2) = inv(inv(K1(x,x2)) + inv(K2(x,x2)) + ... ) Notice that
%	the individual kernels Ki are each allowed to take a different
%	subset of the inputs x, x2 (defined in the field
%	kern.comp{i}.index).
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - the input matrix associated with the rows of the kernel.
%	  X2 - the input matrix associated with the columns of the kernel.
%
%	K = INVCMPNDKERNCOMPUTE(KERN, X) computes the kernel parameters for
%	the kernel that is the inverse of the sum of the precisions of 1 or
%	more kernels given a design matrix of inputs. That is, K(x,x2) =
%	inv(inv(K1(x1,x2)) + inv(K2(x1,x2)) + ... ) Notice that the
%	individual kernels Ki are each allowed to take a different subset of
%	the inputs x, x2 (defined in the field kern.comp{i}.index).
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  X - input data matrix in the form of a design matrix.
%	
%
%	See also
%	INVCMPNDKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, INVCMPNDKERNDIAGCOMPUTE


%	Copyright (c) 2012 Andreas C. Damianou


% !!! TODO: There are several precomputations and tricks to speed this up...!

if nargin > 2
  i = 1;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    k = kernCompute(kern.comp{i}, ...
                         x(:, kern.comp{i}.index), ...
                         x2(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k = kernCompute(kern.comp{i}, x, x2);
  end
  k = pdinv(k); %%%
  for i = 2:length(kern.comp)
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved in the kernel.
      k_cur = kernCompute(kern.comp{i}, ...
                           x(:, kern.comp{i}.index), ...
                           x2(:, kern.comp{i}.index)); 
    else
      % all the data is involved with the kernel.
      k_cur = kernCompute(kern.comp{i}, x, x2);
    end
    k = k + pdinv(k_cur);
  end
else
  i = 1;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    k  = kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k  = kernCompute(kern.comp{i}, x);
  end
  k = pdinv(k); %%%
  for i = 2:length(kern.comp)
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved with the kernel.
      k_cur = kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k_cur = kernCompute(kern.comp{i}, x);
    end
    k = k + pdinv(k_cur);
  end
end

k = pdinv(k);

