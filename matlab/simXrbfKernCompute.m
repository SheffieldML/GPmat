function [K, sK] = simXrbfKernCompute(simKern, rbfKern, t1, t2)

% SIMXRBFKERNCOMPUTE Compute a cross kernel between the SIM and RBF kernels.
% FORMAT
% DESC computes cross kernel terms between SIM and RBF kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
% RETURN sK : normalised matrix (i.e. unscaled version of K not multiplied
% by sqrt(kern.variance) and sqrt(2/simKern.inverseWidth) in the case of
% the un-normalised kernel).
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
% RETURN K : block of values from kernel matrix.
% RETURN sK : normalised matrix (i.e. unscaled version of K not multiplied
% by sqrt(kern.variance) and sqrt(2/simKern.inverseWidth) in the case of
% the un-normalised kernel).
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2008, 2009
%
% MODIFICATIONS : David Luengo, 2009

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
% The SIM kernel implicitly assumes the variance of the RBF kernel
% for f(t) to be 1. To avoid confusion, the same constraint is
% enforced here as well.
if rbfKern.variance ~= 1
  warning(['RBF kernel variance = ' num2str(rbfKern.variance) ...
      '. SIM kernel can only be cross combined with an RBF kernel with variance 1.'])
end
% The normalisation in the SIM kernel arises from the use of the normalised
% version of the RBF kernel. Hence, both have to be normalised or not.
if ~isfield(simKern, 'isNormalised')
    isSimNormalised = false;
else
    isSimNormalised = simKern.isNormalised;
end
if ~isfield(rbfKern, 'isNormalised')
    isRbfNormalised = false;
else
    isRbfNormalised = rbfKern.isNormalised;
end
if isSimNormalised ~= isRbfNormalised
    error('Both the SIM and the RBF kernels have to be either normalised or not.');
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1 = t1 - simKern.delay;
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
sigma = sqrt(2/simKern.inverseWidth);

invSigmaDiffT = 1/sigma*diffT;
halfSigmaD_i = 0.5*sigma*simKern.decay;

%lnPart1 = log(erf(invSigmaDiffT - halfSigmaD_i) + erf(t2Mat/sigma + halfSigmaD_i));
lnPart1 = zeros(size(t1Mat));
warnState = warning('query', 'MATLAB:log:logOfZero');
warning('off', 'MATLAB:log:logOfZero');
if simKern.isStationary == false
    I1 = sign(real(halfSigmaD_i - invSigmaDiffT)) ...
        ~= sign(real(halfSigmaD_i + t2Mat/sigma));
    I2 = (real(halfSigmaD_i - invSigmaDiffT) > 0) & ~I1;
    I3 = ~I2 & ~I1;
    lnPart1(I1) = log(  erfc( real(halfSigmaD_i - invSigmaDiffT(I1))) ...
		   - erfc( real(halfSigmaD_i + t2Mat(I1)/sigma)));
    %lnPart1(I2) = log(  erfc( real(halfSigmaD_i - invSigmaDiffT(I2))) ...
    %		 - erfc( real(halfSigmaD_i + t2Mat(I2)/sigma)));
    %lnPart1(I3) = log(- erfc(-real(halfSigmaD_i - invSigmaDiffT(I3))) ...
    %		 + erfc(-real(halfSigmaD_i + t2Mat(I3)/sigma)));
    lnPart1(I2) = log(  erfcx( real(halfSigmaD_i - invSigmaDiffT(I2))) ...
		   - exp(  real(halfSigmaD_i - invSigmaDiffT(I2)).^2 ...
			   - real(halfSigmaD_i + t2Mat(I2)/sigma).^2) ...
		   .* erfcx( real(halfSigmaD_i + t2Mat(I2)/sigma))) ...
        - real(halfSigmaD_i - invSigmaDiffT(I2)).^2;
    lnPart1(I3) = log(- exp(real(halfSigmaD_i + t2Mat(I3)/sigma).^2 ...
		       - real(halfSigmaD_i - invSigmaDiffT(I3)).^2) ...
            .* erfcx(-real(halfSigmaD_i - invSigmaDiffT(I3))) ...
                + erfcx(-real(halfSigmaD_i + t2Mat(I3)/sigma))) ...
        - real(halfSigmaD_i + t2Mat(I3)/sigma).^2;
    sK = exp(halfSigmaD_i*halfSigmaD_i - simKern.decay*diffT + lnPart1);
else
    lnPart1 = lnDiffErfs(inf, halfSigmaD_i - invSigmaDiffT);
    sK = exp(halfSigmaD_i*halfSigmaD_i - simKern.decay*diffT + lnPart1);
end
warning(warnState.state, 'MATLAB:log:logOfZero');

sK = 0.5 * sK;
if ~isSimNormalised
    sK = sK * sqrt(pi);
    K = sK * sqrt(simKern.variance) * sigma;
else
    K = sK * sqrt(simKern.variance);
end
