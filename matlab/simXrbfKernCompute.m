function K = simXrbfKernCompute(simKern, rbfKern, t1, t2)

% SIMXRBFKERNCOMPUTE Compute a cross kernel between the SIM and RBF kernels.
% FORMAT
% DESC computes cross kernel terms between SIM and RBF kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between SIM and RBF kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN k : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1 = t1 - simKern.delay;
t1Mat = repmat(t1, [1 dim2]);
t2Mat = repmat(t2', [dim1 1]);
diffT = (t1Mat - t2Mat);
sigma = sqrt(2/simKern.inverseWidth);

invSigmaDiffT = 1/sigma*diffT;
halfSigmaD_i = 0.5*sigma*simKern.decay;
K = exp(halfSigmaD_i*halfSigmaD_i)...
    *(exp(-simKern.decay*diffT).*(erf(invSigmaDiffT - halfSigmaD_i) ...
                            + erf(t2Mat/sigma + halfSigmaD_i)));

K = 0.5*sqrt(simKern.variance)*sqrt(rbfKern.variance)*K*sqrt(pi)*sigma;
