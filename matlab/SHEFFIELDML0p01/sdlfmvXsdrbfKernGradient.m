function [g1,g2] = sdlfmvXsdrbfKernGradient(sdlfmvKern, sdrbfKern, t1, ...
    t2, covGrad)

% SDLFMVXSDRBFKERNGRADIENT Gradients cross kernel between a SDLFM and SDRBF
%
%	Description:
%
%	[G1, G2] = SDLFMVXSDRBFKERNGRADIENT(SDLFMVKERN, SDRBFKERN, T1,
%	COVGRAD) computes a cross gradient for a cross kernel between a
%	switching dynamical LFM kernel and a switching dynamical RBF kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	 Arguments:
%	  SDLFMVKERN - the kernel structure associated with the SDLFM kernel
%	   (velocity).
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T1 - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = SDLFMVXSDRBFKERNGRADIENT(SDLFMVKERN, SDRBFKERN, T1, T2,
%	COVGRAD) computes a cross gradient for a cross kernel between a
%	switching dynamical LFM kernel and a switching dynamical RBF kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see lfmKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see lfmKernExtractParam.
%	 Arguments:
%	  SDLFMVKERN - the kernel structure associated with the SDLFM kernel
%	   (velocity).
%	  SDRBFKERN - the kernel structure associated with the SDRBF kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin == 4
    covGrad = t2;
    t2 = t1;
end

[g1, g2] = sdlfmXsdrbfKernGradient(sdlfmvKern, sdrbfKern, t1, t2, covGrad, 'Vel');