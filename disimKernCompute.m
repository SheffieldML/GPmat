function k = disimKernCompute(kern, t, t2)

% DISIMKERNCOMPUTE Compute the DISIM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : disimKernParamInit, kernCompute, kernCreate, disimKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006; Antti Honkela, 2007

% KERN

if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

l = sqrt(2/kern.inverseWidth);

h = disimComputeH(t, t2, kern.di_decay, kern.decay, kern.decay, l);
hp = disimComputeHPrime(t, t2, kern.di_decay, kern.decay, kern.decay, l);
if nargin < 3
  k = h + h' + hp + hp';
else
  h2 = disimComputeH(t2, t, kern.di_decay, kern.decay, kern.decay, l);
  hp2 = disimComputeHPrime(t2, t, kern.di_decay, kern.decay, kern.decay, l);
  k = h + h2' + hp + hp2';
end
k = 0.5*k*sqrt(pi)*l;
k = kern.rbf_variance*kern.di_variance*kern.variance*k;
k = real(k);

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  dim1 = size(t, 1);
  dim2 = size(t2, 1);
  t1Mat = t(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  k = k + kern.initialVariance * kern.variance * ...
      (exp(- kern.di_decay * t1Mat) - exp(- kern.decay * t1Mat)) ./ (kern.decay - kern.di_decay) .* ...
      (exp(- kern.di_decay * t2Mat) - exp(- kern.decay * t2Mat)) ./ (kern.decay - kern.di_decay);
end
