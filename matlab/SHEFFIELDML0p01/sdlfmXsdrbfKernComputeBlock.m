function K = sdlfmXsdrbfKernComputeBlock(lfmKern, rbfKern, t1, t2, ...
    i, j, generalConst)

% SDLFMXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDLFM and SDRBF for i,j
%
%	Description:
%
%	K = SDLFMXSDRBFKERNCOMPUTEBLOCK(LFMKERN, RBFKERN, T1, T2, I, J,
%	GENERALCONSTANT) computes the kernel matrix between a SDLFM kernel
%	function and a SDRBF kernel function in the block specified at
%	indeces i,j. It assumes the computation for a function that describe
%	a position.
%	 Returns:
%	  K - the kernel matrix portion of block i,j
%	 Arguments:
%	  LFMKERN - structure containing parameters for the outputs system
%	  RBFKERN - structure containing parameters for the latent system
%	  T1 - times at which the system 1 is evaluated
%	  T2 - times at which the system 2 is evaluated
%	  I - interval to be evaluated for system 1
%	  J - interval to be evaluated for system 2
%	  GENERALCONSTANT - constants evaluated with
%	   sdlfmKernComputeConstant.m


%	Copyright (c) 2010. Mauricio A. Alvarez


if nargin<7
    j = i;
    generalConst = [];
end

if i==j
    K  = lfmXrbfKernCompute(lfmKern, rbfKern, t1, t2);    
else
    c1 = sdlfmMeanCompute(lfmKern, t1, 'Pos');
    e1 = sdlfmMeanCompute(lfmKern, t1, 'Vel');
    if i>j
        PosRbf = lfmXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
        VelRbf = lfmvXrbfKernCompute(lfmKern, rbfKern, rbfKern.limit, t2);
        if isempty(generalConst{i,j})
            K =  c1*PosRbf + e1*VelRbf;
        else
            K = (generalConst{i,j}(1,1)*c1 + generalConst{i,j}(2,1)*e1)*PosRbf + ...
                (generalConst{i,j}(1,2)*c1 + generalConst{i,j}(2,2)*e1)*VelRbf;
        end
    end
end

