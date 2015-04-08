function [g1, g2] = irbfXrbfKernGradient(irbfKern, rbfKern, t1, t2, covGrad)

% IRBFXRBFKERNGRADIENT Compute gradient between the IRBF and RBF kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between IRBF and RBF kernels for
% the multiple output kernel. 
% ARG irbfKern : the kernel structure associated with the IRBF
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of IRBF kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between IRBF and RBF kernels for
% the multiple output kernel. 
% ARG irbfKern : the kernel structure associated with the IRBF
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of IRBF kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, irbfKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2008, 2009

% KERN

arg{1} = t1;
if nargin < 5
  covGrad = t2;
  t2 = t1;
else
  arg{2} = t2;  
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if irbfKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
% The IRBF kernel implicitly assumes the variance of the RBF kernel
% for f(t) to be 1. To avoid confusion, the same constraint is
% enforced here as well.
if rbfKern.variance ~= 1
  error('IRBF kernel can only be cross combined with an RBF kernel with variance 1.')
end

k = irbfXrbfKernCompute(irbfKern, rbfKern, arg{:});
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
sigma = sqrt(2/irbfKern.inverseWidth);
sigma2 = sigma*sigma;
C_i = sqrt(irbfKern.variance);
D_i = irbfKern.decay;

part2 = exp(-t2Mat.*t2Mat/sigma2-t1Mat*D_i)...
        -exp(-diffT.*diffT/sigma2);

dk_dD = sum(sum((k.*(0.5*sigma2*D_i - diffT) ...
                 + 0.5*C_i*sigma2*part2).*covGrad));

dk_dsigma = sum(sum((k.*(1/sigma+0.5*sigma*D_i*D_i) ...
                     + C_i*sigma* ...
                     ((-diffT/sigma2-D_i/2).*exp(-diffT.*diffT/sigma2)...
                      +(-t2Mat/sigma2+D_i/2).*exp(-t2Mat.*t2Mat/sigma2-t1Mat*D_i))).*covGrad));

dk_dC = sum(sum(k.*covGrad))/C_i;
dk_dRbfVariance = 0;

dk_dinvWidth = -0.5*sqrt(2)/(irbfKern.inverseWidth* ...
                             sqrt(irbfKern.inverseWidth))*dk_dsigma;
dk_dIrbfVariance = dk_dC*0.5/C_i;


% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [dk_dD dk_dinvWidth dk_dIrbfVariance];
g2 = [0 dk_dRbfVariance];


