function [g1, g2] = disimXdisimKernGradient(disimKern1, disimKern2, t1, t2, covGrad)

% DISIMXDISIMKERNGRADIENT Compute a cross gradient between two DISIM kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two disim kernels for the multiple output kernel. 
% ARG disimKern1 : the kernel structure associated with the first DISIM
% kernel.
% ARG disimKern2 : the kernel structure associated with the second DISIM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see disimKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see disimKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between two DISIM kernels for
% the multiple output kernel. 
% ARG disimKern1 : the kernel structure associated with the first DISIM
% kernel.
% ARG disimKern2 : the kernel structure associated with the second DISIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see disimKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see disimKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, disimKernParamInit, disimKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009

% KERN

arg{1}=t1;
if nargin < 5
  covGrad = t2;
  t2 = t1;
else
  arg{2}=t2;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end
if disimKern1.inverseWidth ~= disimKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end
if disimKern1.di_decay ~= disimKern2.di_decay
  error('Kernels cannot be cross combined if they have different driving input decays.');
end
if disimKern1.di_variance ~= disimKern2.di_variance
  error('Kernels cannot be cross combined if they have different driving input variances.');
end
if disimKern1.rbf_variance ~= disimKern2.rbf_variance
  error('Kernels cannot be cross combined if they have different RBF variances.');
end



l = sqrt(2/disimKern1.inverseWidth);
[h1, dh1_ddelta, dh1_dD1, dh1_dD2, dh1_dl] = disimComputeH(t1, t2, disimKern1.di_decay, disimKern1.decay, disimKern2.decay, l);
[hp1, dhp1_ddelta, dhp1_dD1, dhp1_dD2, dhp1_dl] = disimComputeHPrime(t1, t2, disimKern1.di_decay, disimKern1.decay, disimKern2.decay, l);

% Avoid making the expensive call twice with the same arguments
if (all(t1 == t2) && (disimKern1.decay == disimKern2.decay)),
  h2 = h1;
  dh2_ddelta = dh1_ddelta;
  dh2_dD2 = dh1_dD1;
  dh2_dD1 = dh1_dD2;
  dh2_dl = dh1_dl;

  hp2 = hp1;
  dhp2_ddelta = dhp1_ddelta;
  dhp2_dD2 = dhp1_dD1;
  dhp2_dD1 = dhp1_dD2;
  dhp2_dl = dhp1_dl;
else,
  [h2, dh2_ddelta, dh2_dD2, dh2_dD1, dh2_dl] = disimComputeH(t2, t1, disimKern1.di_decay, disimKern2.decay, disimKern1.decay, l);
  [hp2, dhp2_ddelta, dhp2_dD2, dhp2_dD1, dhp2_dl] = disimComputeHPrime(t2, t1, disimKern1.di_decay, disimKern2.decay, disimKern1.decay, l);
end

dK_ddelta = dh1_ddelta + dh2_ddelta' + dhp1_ddelta + dhp2_ddelta';
dK_dD1 = dh1_dD1 + dh2_dD1' + dhp1_dD1 + dhp2_dD1';
dK_dD2 = dh1_dD2 + dh2_dD2' + dhp1_dD2 + dhp2_dD2';
dK_dl = dh1_dl + dh2_dl' + dhp1_dl + dhp2_dl';

C0 = disimKern1.di_variance;
C1 = sqrt(disimKern1.variance);
C2 = sqrt(disimKern2.variance);
C3 = disimKern1.rbf_variance;
K = h1 + h2' + hp1 + hp2';
K = 0.5*K*sqrt(pi);
var2 = C0*C1*C2*C3;
dk_ddelta = (sum(sum(covGrad.*dK_ddelta)))*0.5*sqrt(pi)*l*var2;
dk_dD1 = (sum(sum(covGrad.*dK_dD1)))*0.5*sqrt(pi)*l*var2;
dk_dD2 = (sum(sum(covGrad.*dK_dD2)))*0.5*sqrt(pi)*l*var2;
dk_dl = sum(sum(covGrad.*(dK_dl*0.5*sqrt(pi)*l + K)))*var2;
K = l*K;
dk_dC0 = C1*C2*C3*sum(sum(covGrad.*K));
dk_dC1 = C0*C2*C3*sum(sum(covGrad.*K));
dk_dC2 = C0*C1*C3*sum(sum(covGrad.*K));
dk_dC3 = C0*C1*C2*sum(sum(covGrad.*K));

dk_dDIVariance = dk_dC0;
dk_dDisim1Variance = dk_dC1*0.5/C1;
dk_dDisim2Variance = dk_dC2*0.5/C2;
dk_dRBFVariance = dk_dC3;

dk_dinvWidth = -0.5*sqrt(2)/(disimKern1.inverseWidth* ...
                             sqrt(disimKern1.inverseWidth))*dk_dl;


K = var2*K;

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = [dk_ddelta dk_dinvWidth dk_dDIVariance dk_dD1 dk_dDisim1Variance dk_dRBFVariance];
g2 = [0 0 0 dk_dD2 dk_dDisim2Variance 0];
