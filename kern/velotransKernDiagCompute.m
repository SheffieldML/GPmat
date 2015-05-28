function k = velotransKernDiagCompute(kern, x)

% VELOTRANSKERNDIAGCOMPUTE Compute diagonal of VELOTRANS kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the velocity translate kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix. The last
% column of the input data matrix should be time.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : velotransKernParamInit, kernDiagCompute, kernCreate, velotransKernCompute, translateKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN
  
t = x(:, end);
xPass = x(:, 1:end-1);
xPass = xPass - t*kern.velocity;
k = cmpndKernDiagCompute(kern, xPass);
