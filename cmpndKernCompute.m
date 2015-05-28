function k = cmpndKernCompute(kern, x, x2)

% CMPNDKERNCOMPUTE Compute the CMPND kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the compound
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the compound
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : cmpndKernParamInit, kernCompute, kernCreate, cmpndKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

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
      k  = k + kernCompute(kern.comp{i}, ...
                           x(:, kern.comp{i}.index), ...
                           x2(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k + kernCompute(kern.comp{i}, x, x2);
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
      k  = k + kernCompute(kern.comp{i}, x(:, kern.comp{i}.index));
    else
      % all the data is involved with the kernel.
      k  = k + kernCompute(kern.comp{i}, x);
    end
  end
end
