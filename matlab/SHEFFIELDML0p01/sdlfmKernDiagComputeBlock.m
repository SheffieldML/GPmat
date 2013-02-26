function k = sdlfmKernDiagComputeBlock(lfmKern, t, kyy, kyv, kvy, kvv)

% SDLFMKERNDIAGCOMPUTEBLOCK Diagonal of a SDLFM kernel matrix for block i
%
%	Description:
%
%	K = SDLFMKERNDIAGCOMPUTEBLOCK(LFMKERN, T, KYY, KYV, KVY, KVV)
%	computes the Diagonal of a SDLFM kernel matrix for a particular
%	block It assumes the computation for functions that describe
%	positions (position 1 and position 2).
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



c1 = sdlfmMeanCompute(lfmKern(1), t, 'Pos');
e1 = sdlfmMeanCompute(lfmKern(1), t, 'Vel');

k = kyy*(c1.^2) + (kyv + kvy)*(c1.*e1) + kvv*(e1.^2);

for i=1:length(lfmKern)
    k  = k + lfmKernDiagCompute(lfmKern(i), t);
end