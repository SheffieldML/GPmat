function [g1, g2, g3] = sdlfmXsdrbfKernGradientBlockIEJ(lfmKern, rbfKern, ...
    t1, t2, covGrad, typeGrad)

% SDLFMXSDRBFKERNGRADIENTBLOCKIEJ 
% FORMAT
% DESC computes the gradients of the parameters for output system and
% latent system when i is equal to j, this is, when the kernel function is 
% evaluated at the same switching interval
% ARG lfmKern : structure containing parameters for the output system 
% ARG rbfKern : structure containing parameters for the latent force system
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the corresponding kernel matrix
% ARG typeGrad : specify the mean functions used to compute this part of
% the kernel
% RETURN g1 : gradients of parameters for the system 1
% RETURN g2 : gradients of parameters for the system 2
% RETURN g3 : gradients of switching points
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

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
