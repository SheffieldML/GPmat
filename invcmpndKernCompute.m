function k = invcmpndKernCompute(kern, x, x2)

% INVCMPNDKERNCOMPUTE Compute the INVERSE-PRECISION-CMPND kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the kernel that is the inverse of
% the sum of the precisions of 1 or more kernels, given
% inputs associated with rows and columns. 
% That is, K(x,x2) = inv(inv(K1(x,x2)) + inv(K2(x,x2)) + ... )
% Notice that the individual kernels Ki are each allowed to take a different
% subset of the inputs x, x2 (defined in the field kern.comp{i}.index).
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel parameters for the kernel that is the inverse of
% the sum of the precisions of 1 or more kernels
% given a design matrix of inputs. 
% That is, K(x,x2) = inv(inv(K1(x1,x2)) + inv(K2(x1,x2)) + ... )
% Notice that the individual kernels Ki are each allowed to take a different
% subset of the inputs x, x2 (defined in the field kern.comp{i}.index).
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : invcmpndKernParamInit, kernCompute, kernCreate, invcmpndKernDiagCompute
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN

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

