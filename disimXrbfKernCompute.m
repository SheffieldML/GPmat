function K = disimXrbfKernCompute(disimKern, rbfKern, t1, t2)

% DISIMXRBFKERNCOMPUTE Compute a cross kernel between the DISIM and RBF kernels.
% FORMAT
% DESC computes cross kernel terms between DISIM and RBF kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between DISIM and RBF kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, disimKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009

% KERN

if nargin < 4
  t2 = t1;
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

D_i = disimKern.decay;
delta = disimKern.di_decay;

invLDiffT = 1/l*diffT;
halfLD_i = 0.5*l*D_i;
halfLDelta = 0.5*l*delta;

lnCommon = - log(delta - D_i);
[lnPart1, sign1] = lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l);
[lnPart2, sign2] = lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT);

K = sign1 .* exp(lnCommon + halfLDelta^2 - delta * diffT + lnPart1) ...
    + sign2 .* exp(lnCommon + halfLD_i^2 - D_i * diffT + lnPart2);

K = 0.5*sqrt(disimKern.variance)*sqrt(disimKern.di_variance)*rbfKern.variance*K*sqrt(pi)*l;
K = real(K);
