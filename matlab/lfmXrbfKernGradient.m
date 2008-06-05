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
% MODIFICATIONS : Neil D. Lawrence, 2007; David Luengo, 2008, Mauricio
% Alvarez 2008

% LFM

if nargin < 5
  covGrad = t2;
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern.inverseWidth ~= rbfKern.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

%Kxf = lfmXrbfKernCompute(lfmKern, rbfKern, t1, t2);

m = lfmKern.mass;
D = lfmKern.spring;
C = lfmKern.damper;
S = lfmKern.sensitivity;

alpha = C/(2*m);
omega = sqrt(D/m-alpha^2);

sigma2 = 2/lfmKern.inverseWidth;% Tamporarly changed by Mauricio Alvarez             
sigma = sqrt(sigma2);

gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;

% Initialization of vectors and matrices

T1 = size(t1, 1);
T2 = size(t2, 1);
Tt = repmat(t1, 1, T2);
Tt2 = repmat(transpose(t2), T1, 1);
matGrad = zeros(T1, T2);

% Gradient with respect to m, C and D

for ind = 1:3 % Parameter (m, D or C)
  switch ind
   case 1  % Gradient wrt m
    gradThetaM = 1;
    gradThetaAlpha = -C/(2*(m^2));
    gradThetaOmega = (C^2-2*m*D)/(2*(m^2)*sqrt(4*m*D-C^2));
   case 2  % Gradient wrt D
    gradThetaM = 0;
    gradThetaAlpha = 0;
    gradThetaOmega = 1/sqrt(4*m*D-C^2);
   case 3  % Gradient wrt C
    gradThetaM = 0;
    gradThetaAlpha = 1/(2*m);
    gradThetaOmega = -C/(2*m*sqrt(4*m*D-C^2));
  end
  
  gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
  gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
  
  % Gradient evaluation
  
    if isreal(omega)
        gamma = alpha + j*omega;
        gradThetaGamma = gradThetaAlpha + j*gradThetaOmega;
        matGrad(:,:) = -(sigma*sqrt(pi)*S/(2*m*omega)) ...
            * imag(lfmGradientUpsilon(gamma,sigma2,gradThetaGamma,Tt,Tt2) ...
            - (gradThetaM/m + gradThetaOmega/omega) ...
            * lfmComputeUpsilon(gamma,sigma2,Tt,Tt2));
    else
        gamma1 = alpha + j*omega;
        gamma2 = alpha - j*omega;
        gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
        gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
        matGrad(:,:) = (sigma*sqrt(pi)*S/(j*4*m*omega)) ...
            * (lfmGradientUpsilon(gamma2,sigma2,gradThetaGamma2,Tt,Tt2) ...
            - lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma1,Tt,Tt2) ...
            - (gradThetaM/lfmKern.mass + gradThetaOmega/omega) ...
            * (lfmComputeUpsilon(gamma2,sigma2,Tt,Tt2) ...
            - lfmComputeUpsilon(gamma1,sigma2,Tt,Tt2)));
    end

    g1(ind) = sum(sum(matGrad.*covGrad));

end

% Gradient with respect to sigma

if isreal(omega)
    gamma = alpha + j*omega;
    matGrad(:,:) = -(sqrt(pi)*S/(2*m*omega)) ...
        * imag(lfmComputeUpsilon(gamma,sigma2,Tt,Tt2) ...
        + sigma*lfmGradientSigmaUpsilon(gamma,sigma2,Tt,Tt2));
else
    gamma1 = alpha + j*omega;
    gamma2 = alpha - j*omega;
    matGrad(:,:) = (sqrt(pi)*S/(j*4*m*omega)) ...
        *(lfmComputeUpsilon(gamma2,sigma2,Tt,Tt2) ...
        - lfmComputeUpsilon(gamma1,sigma2,Tt,Tt2) ...
        + sigma*(lfmGradientSigmaUpsilon(gamma2,sigma2,Tt,Tt2) ...
        - lfmGradientSigmaUpsilon(gamma1,sigma2,Tt,Tt2)));
end;

% g1(4) = sum(sum(matGrad.*Kxf.*covGrad))*(-(sigma^3)/4); % temporarly introduced by Mauricio Alvarez
g1(4) = sum(sum(matGrad.*covGrad))*(-(sigma^3)/4); % temporarly introduced by Mauricio Alvarez
g2(1) = g1(4);

% Gradient with respect to S

if isreal(omega)
    gamma = alpha + j*omega;
    matGrad(:,:) = -(sqrt(pi)*sigma/(2*m*omega)) ...
        * imag(lfmComputeUpsilon(gamma,sigma2,Tt,Tt2));
else
    gamma1 = alpha + j*omega;
    gamma2 = alpha - j*omega;
    matGrad(:,:) = (sqrt(pi)*sigma/(j*4*m*omega)) ...
        *(lfmComputeUpsilon(gamma2,sigma2,Tt,Tt2) - ...
          lfmComputeUpsilon(gamma1,sigma2,Tt,Tt2));
end;

g1(5) = sum(sum(matGrad.*covGrad));

% Gradient with respect to the "variance" of the RBF
g2(1) = 0; % Otherwise is counted twice, temporarly changed by Mauricio Alvarez
g2(2) = 0;
