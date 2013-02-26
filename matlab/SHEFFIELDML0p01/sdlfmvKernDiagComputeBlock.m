function k = sdlfmvKernDiagComputeBlock(lfmKern, t, kyy, kyv, kvy, kvv)

% SDLFMVKERNDIAGCOMPUTEBLOCK Diagonal of a SDLFM kernel matrix for block i
%
%	Description:
%
%	K = SDLFMVKERNDIAGCOMPUTEBLOCK(LFMKERN, T, KYY, KYV, KVY, KVV)
%	computes the diagonal of a SDLFM kernel matrix for a particular
%	block It assumes the computation for functions that describe
%	velocities
%	 Returns:
%	  K - the diagonal of the kernel matrix portion of a block
%	 Arguments:
%	  LFMKERN - structure containing parameters for the system
%	  T - times at which the system is evaluated
%	  KYY - covariance for the initial conditions between position 1 and
%	   position 2 at block i,j
%	  KYV - covariance for the initial conditions between position 1 and
%	   velocity 2 at block i,j
%	  KVY - covariance for the initial conditions between velocity 1 and
%	   position 2 at block i,j
%	  KVV - covariance for the initial conditions between velocity 1 and
%	   velocity 2 at block i,j


%	Copyright (c) 2010. Mauricio A. Alvarez


g1 = sdlfmvMeanCompute(lfmKern(1), t, 'Pos');
h1 = sdlfmvMeanCompute(lfmKern(1), t, 'Vel');

k = kyy*(g1.^2) + (kyv + kvy)*(g1.*h1) + kvv*(h1.^2);

for i=1:length(lfmKern)
    k  = k + lfmvKernDiagCompute(lfmKern(i), t);
end