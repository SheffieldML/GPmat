function [g1, g2] = lfmXlfmKernGradient(lfmKern1, lfmKern2, t, t2, covGrad)

% LFMXLFMKERNGRADIENT Compute a cross gradient between two LFM kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two lfm kernels for the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between two LFM kernels for
% the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmKernParamInit, lfmKernExtractParam
%
% COPYRIGHT : David Luengo, 2007
% 
% MODIFICATIONS : Neil D. Lawrence, 2007

% KERN

  
if nargin < 5
  covGrad = t2;
  t2 = t;
end
if size(t, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

Kxx = lfmXlfmKernCompute(lfmKern1, lfmKern2, t, t2);

% Get length scale out.
sigma2 = 2/lfmKern1.inverseWidth;
sigma = sqrt(sigma2);


if false
    autoCorr = true;
    if ind_par
        ind_par = 0;
        warning(['Inconsistency in the input parameters. ' ...
                 'Assuming that Kxx is an autocorrelation matrix.']);
    end
else
    autoCorr = false;
end

C = [lfmKern1.damper lfmKern2.damper];
D = [lfmKern1.spring lfmKern2.spring];
m = [lfmKern1.mass lfmKern2.mass];

% Parameters of the simulation
alpha = C./(2*m);
omega = sqrt(D./m-alpha.^2);

sigma = sqrt(sigma2);

gamma1 = alpha + j*omega;
gamma2 = alpha - j*omega;

% Choosing the right gradients for m, omega, gamma1 and gamma2

for ind_theta = 1:3
  switch ind_theta
   case 1 
    gradThetaM = 1;
    gradThetaAlpha = -C./(2*(m.^2));
    gradThetaOmega = (C.^2-2*m.*D)./(2*(m.^2).*sqrt(4*m.*D-C.^2));
   case 2
    gradThetaM = 0;
    gradThetaAlpha = 1./(2*m);
    gradThetaOmega = -C./(2*m.*sqrt(4*m.*D-C.^2));
   case 3
    gradThetaM = 0;
    gradThetaAlpha = 0;
    gradThetaOmega = 1./sqrt(4*m.*D-C.^2);
  end
  
  gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
  gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
  
  for ind_par = 0:1
    if autoCorr
      matInd = ones(8,2);
    else
      matInd = repmat([ind_par 1+ind_par; 1+ind_par ind_par],4,1);
    end
    
    gradThetaGamma = [gradThetaGamma2(2) gradThetaGamma1(1); ...
                      gradThetaGamma1(1) gradThetaGamma2(2); ...
                      gradThetaGamma1(2) gradThetaGamma2(1); ...
                      gradThetaGamma2(1) gradThetaGamma1(2); ...
                      gradThetaGamma2(2) gradThetaGamma2(1); ...
                      gradThetaGamma2(1) gradThetaGamma2(2); ...
                      gradThetaGamma1(2) gradThetaGamma1(1); ...
                      gradThetaGamma1(1) gradThetaGamma1(2)].*matInd;
    
    % Gradient evaluation
    matGrad = (sigma*lfmKern1.sensitivity*lfmKern2.sensitivity*sqrt(pi) ...
               /(8*(lfmKern1.mass*lfmKern2.mass*prod(omega))) ...
               *(lfmGradientH(gamma2(2),gamma1(1),sigma2, ...
                              gradThetaGamma(1,:),t,t2) ...
                 + transpose(lfmGradientH(gamma1(1),gamma2(2),sigma2, ...
                                          gradThetaGamma(2,:),t2,t)) ...
                 + lfmGradientH(gamma1(2),gamma2(1),sigma2, ...
                                gradThetaGamma(3,:),t,t2) ...
                 + transpose(lfmGradientH(gamma2(1),gamma1(2),sigma2, ...
                                          gradThetaGamma(4,:),t2,t)) ...
                 - lfmGradientH(gamma2(2),gamma2(1),sigma2, ...
                                gradThetaGamma(5,:),t,t2) ...
                 - transpose(lfmGradientH(gamma2(1),gamma2(2),sigma2, ...
                                          gradThetaGamma(6,:),t2,t)) ...
                 - lfmGradientH(gamma1(2),gamma1(1),sigma2, ...
                                gradThetaGamma(7,:),t,t2) ...
                 - transpose(lfmGradientH(gamma1(1),gamma1(2),sigma2, ...
                                          gradThetaGamma(8,:),t2,t))) ...
               - (1+autoCorr)*Kxx*(gradThetaM/m(1+ind_par) ...
                                   + gradThetaOmega(1+ind_par)/omega(1+ind_par)));
    
    if ind_par == 0
      g1(ind_theta) = sum(sum(matGrad.*covGrad));
    else
      g2(ind_theta) = sum(sum(matGrad.*covGrad));
    end
  end
end

matGrad = Kxx/sigma ...
          + (sigma*lfmKern1.sensitivity*lfmKern2.sensitivity*sqrt(pi)/(8*prod(m)*prod(omega)))*...
          (lfmGradientSigmaH(gamma2(1),gamma1(2),sigma2,t,t2) + ...
           transpose(lfmGradientSigmaH(gamma1(2),gamma2(1),sigma2,t2,t)) + ...
           lfmGradientSigmaH(gamma1(1),gamma2(2),sigma2,t,t2) + ...
           transpose(lfmGradientSigmaH(gamma2(2),gamma1(1),sigma2,t2,t)) - ...
           lfmGradientSigmaH(gamma2(1),gamma2(2),sigma2,t,t2) - ...
           transpose(lfmGradientSigmaH(gamma2(2),gamma2(1),sigma2,t2,t)) - ...
           lfmGradientSigmaH(gamma1(1),gamma1(2),sigma2,t,t2) - ...
           transpose(lfmGradientSigmaH(gamma1(2),gamma1(1),sigma2,t2,t)));

g1(4) = sum(sum(matGrad.*Kxx));
g1(4) = -g1(4)/4*(sigma*sigma*sigma);
g2(4) = g1(4);

lfmKern1.inverseWidth= 2/sigma*sigma;
g1(5) = sum(sum(Kxx/lfmKern1.sensitivity.*covGrad));
g2(5) = sum(sum(Kxx/lfmKern2.sensitivity.*covGrad));


