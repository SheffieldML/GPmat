function [g1, g2, g3] = sdlfmaXsdrbfKernGradientBlock(lfmKern, rbfKern, t1, t2, ...
    i, j, generalConst, generalConstGrad, covGrad)

% SDLFMAXSDRBFKERNGRADIENTBLOCK Gradients of the parameters in block i,j
%
%	Description:
%
%	[G1, G2, G3] = SDLFMAXSDRBFKERNGRADIENTBLOCK(LFMKERN, RBFKERN, T1,
%	T2, I, J, GENERALCONSTANT, GENERALCONSTGRAD, COVGRAD) computes the
%	kernel parameters gradients for the SDLFMXSDRBF kernel function in
%	the block specified at indeces i,j. Assumes the output system is an
%	acceleration.
%	 Returns:
%	  G1 - gradients of parameters for the system 1
%	  G2 - gradients of parameters for the system 2
%	  G3 - gradients of switching points
%	 Arguments:
%	  LFMKERN - structure containing parameters for the output system
%	  RBFKERN - structure containing parameters for the latent system
%	  T1 - times at which the system 1 is evaluated
%	  T2 - times at which the system 2 is evaluated
%	  I - interval to be evaluated for system 1
%	  J - interval to be evaluated for system 2
%	  GENERALCONSTANT - constants evaluated with
%	   sdlfmKernComputeConstant.m
%	  GENERALCONSTGRAD - derivatives of the constants computed with
%	   sdlfmKernGradientConstant.m
%	  COVGRAD - partial derivatives of the objective function wrt
%	   portion of the kernel matrix in block i,j


%	Copyright (c) 2010. Mauricio A. Alvarez


if nargin<7
    j = i;
    generalConst = [];
end

if i==j
    [g1, g2, g3t] = sdlfmXsdrbfKernGradientBlockIEJ(lfmKern, rbfKern, t1, t2, ...
        covGrad, {'lfma', 'lfmj'});
    g3(i) = g3t;
else
    [g1, g2, g3] = sdlfmXsdrbfKernGradientBlockIGJ(lfmKern, rbfKern, ...
        t1, t2, i, j, generalConst, generalConstGrad, covGrad, {'lfma', 'lfmj'});
end