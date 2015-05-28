function [g1, g2] = simXrbfKernGradient(simKern, rbfKern, t1, t2, covGrad)

% SIMXRBFKERNGRADIENT Compute gradient between the SIM and RBF kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between SIM and RBF kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SIM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between SIM and RBF kernels for
% the multiple output kernel. 
% ARG simKern : the kernel structure associated with the SIM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of SIM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, simKernParamInit, rbfKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2008, 2009
%
% MODIFICATIONS : David Luengo, 2009

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

[k, sK] = simXrbfKernCompute(simKern, rbfKern, arg{:});
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1 = t1 - simKern.delay;
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
sigma = sqrt(2/simKern.inverseWidth);
sigma2 = sigma*sigma;

if isfield(simKern, 'isNegativeS') && (simKern.isNegativeS == true)
    C_i = simKern.sensitivity;
else
    C_i = sqrt(simKern.variance);
end

D_i = simKern.decay;
N_i = 1/sqrt(2*pi*simKern.inverseWidth);

if (simKern.isStationary == false)
    part2 = exp(-t2Mat.*t2Mat/sigma2-t1Mat*D_i) - exp(-diffT.*diffT/sigma2);
    if ~isSimNormalised
        dk_dD = sum(sum((k.*(0.5*sigma2*D_i - diffT) + 0.5*C_i*sigma2*part2).*covGrad));
        dk_dsigma = sum(sum((C_i * (1+0.5*sigma2*D_i*D_i) * sK ...
            + C_i * ((-diffT/sigma-0.5*sigma*D_i).*exp(-diffT.*diffT/sigma2) ...
                + (-t2Mat/sigma+0.5*sigma*D_i).*exp(-t2Mat.*t2Mat/sigma2-t1Mat*D_i))).*covGrad));
    else
        dk_dD = sum(sum((k.*(0.5*sigma2*D_i - diffT) + 0.5*C_i*sigma/sqrt(pi)*part2).*covGrad));
        dk_dsigma = sum(sum((0.5 * sigma * C_i * D_i * D_i * sK ...
            + C_i/sqrt(pi) * ((-diffT/sigma2-0.5*D_i).*exp(-diffT.*diffT/sigma2) ...
                +(-t2Mat/sigma2+0.5*D_i).*exp(-t2Mat.*t2Mat/sigma2-t1Mat*D_i))).*covGrad));
    end
else
    if ~isSimNormalised
        dk_dD = sum(sum((k .* (0.5*sigma2*D_i - diffT) ...
            - 0.5 * C_i * sigma2 * exp(-diffT.*diffT/sigma2)) .* covGrad));
        dk_dsigma = sum(sum((C_i * (1+0.5*sigma2*D_i*D_i) * sK...
            - C_i * (diffT/sigma+0.5*sigma*D_i) .* exp(-diffT.*diffT/sigma2)) .* covGrad));
    else
        dk_dD = sum(sum((k .* (0.5*sigma2*D_i - diffT) ...
            - 0.5 * C_i * sigma/sqrt(pi) * exp(-diffT.*diffT/sigma2)) .* covGrad));
        dk_dsigma = sum(sum((0.5 * sigma * C_i * D_i * D_i * sK...
            - C_i/sqrt(pi) * (diffT/sigma2+0.5*D_i) .* exp(-diffT.*diffT/sigma2)) .* covGrad));
    end
end

if ~isSimNormalised
    dk_dC = sigma * sum(sum(sK.*covGrad));
else
    dk_dC = sum(sum(sK.*covGrad));
end

dk_dRbfVariance = 0;

dk_dinvWidth = -0.5*sqrt(2)/(simKern.inverseWidth* ...
                             sqrt(simKern.inverseWidth))*dk_dsigma;

if isfield(simKern, 'isNegativeS') && (simKern.isNegativeS == true)
    dk_dSimVariance = dk_dC;
else
    dk_dSimVariance = dk_dC*0.5/C_i;
end
                         
% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = real([dk_dD dk_dinvWidth dk_dSimVariance]);
g2 = [0 dk_dRbfVariance];
