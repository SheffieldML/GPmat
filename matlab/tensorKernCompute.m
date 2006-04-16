function k = tensorKernCompute(kern, x, x2)

% TENSORKERNCOMPUTE Compute the tensor kernel given the parameters and X.
% FORMAT
% DESC computes a kernel matrix for a tensor kernel type given an
% input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : computed elements of the kernel structure.
%
% FORMAT 
% DESC computes a kernel matrix for the tensor kernel type given
% two input data matrices, one for the rows and one for the columns.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the rows of the kernel).
% ARG X2 : second input matrix to the kernel computation (forms the columns of the kernel).
% RETURN K : computed elements of the kernel structure.
%
% SEEALSO : kernCompute, tensorKernCreate, tensorKernDiagCompute

% KERN

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
  for i = 2:length(kern.comp)
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved in the kernel.
      k  = k.*kernCompute(kern.comp{i}, ...
                          x(:, kern.comp{i}.index), ...
                          x2(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k.*kernCompute(kern.comp{i}, x, x2);
    end
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
  for i = 2:length(kern.comp)
    if ~isempty(kern.comp{i}.index)
      % only part of the data is involved with the kernel.
      k  = k.*kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k.*kernCompute(kern.comp{i}, x);
    end
  end
end
