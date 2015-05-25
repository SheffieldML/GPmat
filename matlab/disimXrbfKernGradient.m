function [g1, g2] = disimXrbfKernGradient(disimKern, rbfKern, t1, t2, covGrad)

% DISIMXRBFKERNGRADIENT Compute gradient between the DISIM and RBF kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between DISIM and RBF kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of DISIM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between DISIM and RBF kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of DISIM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, disimKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009

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
if disimKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
if disimKern.rbf_variance ~= rbfKern.variance
  error('Kernels cannot be cross combined if they have different RBF variances.');
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
l = sqrt(2/disimKern.inverseWidth);
l2 = l*l;
C_0 = sqrt(disimKern.di_variance);
C_i = sqrt(disimKern.variance);
C_j = rbfKern.variance;
D_i = disimKern.decay;
delta = disimKern.di_decay;

invLDiffT = 1/l*diffT;
halfLD_i = 0.5*l*D_i;
halfLDelta = 0.5*l*delta;

prefact = C_0 * C_i * C_j * sqrt(pi)/2 * l;

lnCommon = - log(delta - D_i);
[lnPart1, sign1] = lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l);
[lnPart2, sign2] = lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT);

lnFact1 = halfLDelta^2 - delta * diffT;
lnFact2 = halfLD_i^2 - D_i * diffT;


k = sign1 .* exp(lnCommon + lnFact1 + lnPart1) ...
    + sign2 .* exp(lnCommon + lnFact2 + lnPart2);

k = 0.5*sqrt(disimKern.variance)*sqrt(disimKern.di_variance)*rbfKern.variance*k*sqrt(pi)*l;



[dlnPart1, m1] = gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, ...
				l/2, l/2);
[dlnPart2, m2] = gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, ...
				l/2, l/2);

dK_dD = k .* (1./(delta-D_i)) ...
	+ prefact ...
	.* ((l*halfLD_i - diffT) .* sign2 .* exp(lnCommon + lnFact2 + lnPart2) ...
	    + dlnPart2 .* exp(lnCommon + lnFact2 - m2));
dk_dD = sum(sum(dK_dD.*covGrad));

dK_ddelta = k .* (-1./(delta-D_i)) ...
    + prefact ...
    .* ((l*halfLDelta - diffT) .* sign1 .* exp(lnCommon + lnFact1 + lnPart1) ...
	+ dlnPart1 .* exp(lnCommon + lnFact1 - m1));
dk_ddelta = sum(sum(dK_ddelta.*covGrad));

[dlnPart1, m1] = gradLnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l, ...
				delta/2 + invLDiffT/l, delta/2 - t2Mat/l2);
[dlnPart2, m2] = gradLnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT, ...
				D_i/2 - t2Mat/l2, D_i/2 + invLDiffT/l);

dK_dl = k/l ...
	+ prefact ... 
	.* (delta*halfLDelta .* sign1 .* exp(lnCommon + lnFact1 + lnPart1) ...
	    + D_i*halfLD_i .* sign2 .* exp(lnCommon + lnFact2 + lnPart2) ...
	    + dlnPart1 .* exp(lnCommon + lnFact1 - m1) ...
	    + dlnPart2 .* exp(lnCommon + lnFact2 - m2));
dk_dl = sum(sum(dK_dl.*covGrad));

dk_dC_i = sum(sum(k.*covGrad))/C_i;
dk_dC_0 = sum(sum(k.*covGrad))/C_0;
dk_dRbfVariance = sum(sum(k.*covGrad))/rbfKern.variance;

dk_dinvWidth = -0.5*sqrt(2)/(disimKern.inverseWidth* ...
                             sqrt(disimKern.inverseWidth))*dk_dl;
dk_dDisimVariance = dk_dC_i*0.5/C_i;
dk_dDIVariance = dk_dC_0*0.5/C_0;


% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = real([dk_ddelta dk_dinvWidth dk_dDIVariance dk_dD dk_dDisimVariance 0]);
g2 = real([0 dk_dRbfVariance]);


