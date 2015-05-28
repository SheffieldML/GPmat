function K = disimXsimKernCompute(disimKern, simKern, t1, t2)

% DISIMXSIMKERNCOMPUTE Compute a cross kernel between DISIM and SIM kernels.
% FORMAT
% DESC computes cross kernel terms between DISIM and SIM kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between DISIM and SIM kernels for
% the multiple output kernel. 
% ARG disimKern : the kernel structure associated with the DISIM
% kernel.
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, disimKernParamInit, simKernParamInit
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
if disimKern.inverseWidth ~= simKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
if disimKern.di_decay ~= simKern.decay
  error('Kernels cannot be cross combined if they have different driving input decays.');
end
if disimKern.di_variance ~= simKern.variance
  error('Kernels cannot be cross combined if they have different driving input variances.');
end
if disimKern.rbf_variance ~= 1
  warning('KERN:simRBFVariance', ...
	  'Warning: copying non-unit RBF variance from DISIM to SIM kernel')
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

lnCommon1 = - log(2*delta) + halfLDelta.^2;

lnFact1 = log(2 * delta) - log(delta^2 - D_i^2) -delta * t2Mat - D_i * t1Mat;
lnPart1 = lnDiffErfs(halfLDelta - t2Mat/l, halfLDelta);

lnFact2 = -delta * (t1Mat+t2Mat) - log(delta - D_i);
lnPart2a = lnDiffErfs(halfLDelta, halfLDelta - t1Mat/l);
lnPart2b = lnDiffErfs(halfLDelta, halfLDelta - t2Mat/l);

lnFact3 = delta * diffT - log(delta + D_i);
lnPart3 = lnDiffErfs(halfLDelta + t1Mat/l, halfLDelta + invLDiffT);

lnFact4 = -delta*diffT - log(delta - D_i);
lnPart4 = lnDiffErfs(halfLDelta - invLDiffT, halfLDelta + t2Mat/l);

lnCommon2 = - log(delta^2 - D_i^2) - delta * t2Mat - D_i * t1Mat + halfLD_i^2;
lnPart5 = lnDiffErfs(halfLD_i - t1Mat/l, halfLD_i);

lnFact6 = (D_i + delta) * t2Mat;
lnPart6 = lnDiffErfs(halfLD_i + t2Mat/l, halfLD_i - invLDiffT);

K = exp(lnCommon1 + lnFact1 + lnPart1) ...
    +exp(lnCommon1 + lnFact2 + lnPart2a) ...
    +exp(lnCommon1 + lnFact2 + lnPart2b) ...
    +exp(lnCommon1 + lnFact3 + lnPart3) ...
    +exp(lnCommon1 + lnFact4 + lnPart4) ...
    +exp(lnCommon2           + lnPart5) ...
    +exp(lnCommon2 + lnFact6 + lnPart6);
K = 0.5*K*sqrt(pi)*l;
K = disimKern.rbf_variance*disimKern.di_variance*sqrt(disimKern.variance)*K;
K = real(K);

if isfield(disimKern, 'gaussianInitial') && disimKern.gaussianInitial && ...
  isfield(simKern, 'gaussianInitial') && simKern.gaussianInitial,
  if disimKern.initialVariance ~= simKern.initialVariance
    error('Kernels cannot be cross combined if they have different initial variances.');
  end
  
  dim1 = size(t1, 1);
  dim2 = size(t2, 1);
  t1Mat = t1(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  delta = disimKern.di_decay;
  D = disimKern.decay;
  
  K = K + disimKern.initialVariance * sqrt(disimKern.variance) * ...
      (exp(-delta * t1Mat) - exp(-D * t1Mat)) ./ (D - delta) .* ...
      exp(-delta * t2Mat);
end
