function k = sdlfmKernDiagComputeBlock(lfmKern, t, kyy, kyv, kvy, kvv)

% SDLFMKERNDIAGCOMPUTEBLOCK Diagonal of a SDLFM kernel matrix for block i
% FORMAT
% DESC computes the Diagonal of a SDLFM kernel matrix for a particular block  
% It assumes the computation for functions that describe positions 
% (position 1 and position 2).
% ARG lfmKern : structure containing parameters for the system
% ARG t : times at which the system is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% ARG kyv : covariance for the initial conditions between position 1 and
% velocity 2 at block i,j
% ARG kvy : covariance for the initial conditions between velocity 1 and
% position 2 at block i,j
% ARG kvv : covariance for the initial conditions between velocity 1 and
% velocity 2 at block i,j
% RETURN k : the diagonal of the kernel matrix portion of a block
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN


c1 = sdlfmMeanCompute(lfmKern(1), t, 'Pos');
e1 = sdlfmMeanCompute(lfmKern(1), t, 'Vel');

k = kyy*(c1.^2) + (kyv + kvy)*(c1.*e1) + kvv*(e1.^2);

for i=1:length(lfmKern)
    k  = k + lfmKernDiagCompute(lfmKern(i), t);
end
