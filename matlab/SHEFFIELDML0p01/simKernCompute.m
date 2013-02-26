function [k, sk] = simKernCompute(kern, t, t2)

% SIMKERNCOMPUTE Compute the SIM kernel given the parameters and X.
%
%	Description:
%
%	[K, SK] = SIMKERNCOMPUTE(KERN, T1, T2) computes the kernel
%	parameters for the single input motif kernel given inputs associated
%	with rows and columns.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - unscaled kernel matrix (i.e. only 0.5 times the sum of h's
%	   part).
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T1 - the input matrix associated with the rows of the kernel.
%	  T2 - the input matrix associated with the columns of the kernel.
%
%	[K, SK] = SIMKERNCOMPUTE(KERN, T) computes the kernel matrix for the
%	single input motif kernel given a design matrix of inputs.
%	 Returns:
%	  K - the kernel matrix computed at the given points.
%	  SK - unscaled kernel matrix (i.e. only 0.5 times the sum of h's
%	   part).
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.
%	
%	
%	
%	
%
%	See also
%	SIMKERNPARAMINIT, KERNCOMPUTE, KERNCREATE, SIMKERNDIAGCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence


%	With modifications by David Luengo 2009


%	With modifications by Mauricio Alvarez 2009


%	With modifications by Antti Honkela 2009


if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);

if (kern.isStationary == false)
    h = simComputeH(t, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
else
    h = simComputeHStat(t, t2, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
end
if nargin < 3
  sk = 0.5 * (h + h');
else
  if (kern.isStationary == false)
      h2 = simComputeH(t2, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
  else
      h2 = simComputeHStat(t2, t, kern.decay, kern.decay, kern.delay, kern.delay, sigma);
  end
  sk = 0.5 * (h + h2');
end
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    k = sk;
else
    k = sqrt(pi)*sigma*sk;
end

if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    k = (kern.sensitivity*kern.sensitivity)*k;
else
    k = kern.variance*k;
end

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  dim1 = size(t, 1);
  dim2 = size(t2, 1);
  t1Mat = t(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  k = k + kern.initialVariance * exp(- kern.decay * (t1Mat + t2Mat));
end
