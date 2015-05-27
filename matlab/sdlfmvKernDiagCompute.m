function k = sdlfmvKernDiagCompute(sdlfmvKern, t, covIC)

% SDLFMVKERNDIAGCOMPUTE Compute diagonal of a SDLFMV kernel.
% FORMAT
% DESC computes the diagonal of the velocities of the kernel matrix for 
% the switching dynamical latent force model kernel given a design matrix 
% of inputs.
% ARG sdlfmvKern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% ARG covIC : covariance of the initial position
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sdlfmvKernParamInit, kernDiagCompute
%
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

k = sdlfmKernDiagCompute(sdlfmvKern, t, covIC, 'VelVel');
