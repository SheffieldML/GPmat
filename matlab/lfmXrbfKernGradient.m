function [g1, g2] = lfmXrbfKernGradient(lfmKern, rbfKern, t1, t2, covGrad)

% LFMXRBFKERNGRADIENT Compute gradient between the LFM and RBF kernels.
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between LFM and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of LFM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% FORMAT
% DESC computes the gradient of an objective function with respect
% to cross kernel terms between LFM and RBF kernels for
% the multiple output kernel. 
% ARG lfmKern : the kernel structure associated with the LFM
% kernel.
% ARG rbfKern : the kernel structure associated with the RBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN g1 : gradient of objective function with respect to kernel
% parameters of LFM kernel.
% RETURN g2 : gradient of objective function with respect to kernel
% parameters of RBF kernel.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmKernParamInit, rbfKernParamInit
%
% COPYRIGHT : David Luengo, 2007
%
% MODIFICATIONS : Neil D. Lawrence, 2007

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
if lfmKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

Kxf = lfmXrbfKernCompute(lfmKern, rbfKern, arg{:});

C = lfmKern.damper;
D = lfmKern.spring;
m = lfmKern1.mass;

alpha = C/(2*m);
omega = sqrt(D/m-alpha^2);

sigma = sqrt(sigma2);

gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;

% Initialization of vectors and matrices

matGrad = zeros(T,T2);

% Choosing the right gradients for m, omega, gamma1 and gamma2

ParList = upper(strvcat('m', 'C', 'D'));
ind = strmatch(upper(par),ParList,'exact');
if isempty(ind)
    error('Gradient parameter not valid');
end;

for ind = 1:3
  switch ind
   case 1
    gradThetaM = 1;
    gradThetaAlpha = -C/(2*(m^2));
    gradThetaOmega = (C^2-2*m*D)/(2*(m^2)*sqrt(4*m*D-C^2));
   case 2
    gradThetaM = 0;
    gradThetaAlpha = 1/(2*m);
    gradThetaOmega = -C/(2*m*sqrt(4*m*D-C^2));
   case 3
    gradThetaM = 0;
    gradThetaAlpha = 0;
    gradThetaOmega = 1/sqrt(4*m*D-C^2);
  end
  
  gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
  gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
  
  % Gradient evaluation
  
  matGrad = -Kxf*(gradThetaM/m + gradThetaOmega/omega) + ...
            (sigma*S*sqrt(pi)/(j*4*m*omega))*...
            (lfmGradientUpsilon(gamma2,sigma2,gradThetaGamma2,t,t2) - ...
             lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma1,t,t2));
  
  g1(ind) = sum(sum(matGrad.*covGrad));

end
g1(4) = sum(sum(covGrad.*Kxf/lfmKern.sensitivity));
