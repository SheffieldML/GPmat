function [g1, g2, g3] = sdlfmXsdrbfKernGradientBlockIEJ(lfmKern, rbfKern, ...
    t1, t2, covGrad, typeGrad)

% SDLFMXSDRBFKERNGRADIENTBLOCKIEJ
%
%	Description:
%
%	[G1, G2, G3] = SDLFMXSDRBFKERNGRADIENTBLOCKIEJ(LFMKERN, RBFKERN, T1,
%	T2, COVGRAD, TYPEGRAD) computes the gradients of the parameters for
%	output system and latent system when i is equal to j, this is, when
%	the kernel function is evaluated at the same switching interval
%	 Returns:
%	  G1 - gradients of parameters for the system 1
%	  G2 - gradients of parameters for the system 2
%	  G3 - gradients of switching points
%	 Arguments:
%	  LFMKERN - structure containing parameters for the output system
%	  RBFKERN - structure containing parameters for the latent force
%	   system
%	  T1 - times at which the system 1 is evaluated
%	  T2 - times at which the system 2 is evaluated
%	  COVGRAD - partial derivatives of the objective function wrt
%	   portion of the corresponding kernel matrix
%	  TYPEGRAD - specify the mean functions used to compute this part of
%	   the kernel


%	Copyright (c) 2010. Mauricio A. Alvarez


if nargin < 6
    typeGrad  = {'lfm', 'lfmv'};
end

fhandleParam = str2func([typeGrad{1} 'XrbfKernGradient']); 
fhandleSwitching1 = str2func([typeGrad{2} 'XrbfKernCompute']);
fhandleSwitching2 = str2func([typeGrad{1} 'XrbfvKernCompute']);

% Get gradients wrt the parameters

[g1, g2] = fhandleParam(lfmKern, rbfKern, t1, t2, covGrad);

gspKern1 = fhandleSwitching1(lfmKern, rbfKern, t1, t2);
gspKern2 = fhandleSwitching2(lfmKern, rbfKern, t1, t2);

g3 = -sum(sum((gspKern1 + gspKern2).*covGrad));