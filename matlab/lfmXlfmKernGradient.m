function [g1, g2] = lfmXlfmKernGradient(lfmKern1, lfmKern2, t1, t2, covGrad)

% LFMXLFMKERNGRADIENT Compute a cross gradient between two LFM kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two lfm kernels for the multiple output kernel. 
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
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
% COPYRIGHT : David Luengo, 2007, 2008
% 
% MODIFICATIONS : Neil D. Lawrence, 2007
%
% MODIFICATIONS : Mauricio Alvarez, 2008

% KERN

if nargin < 5
  covGrad = t2;
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

%Kxx = lfmXlfmKernCompute(lfmKern1, lfmKern2, t1, t2);

% Check whether lfmKern1 is the same as lfmKern2 (autocorrelation) or not.
% For the time being we consider that both kernels are equal if all their
% parameters are the same

par1 = kernExtractParam(lfmKern1);
par2 = kernExtractParam(lfmKern2);
if (par1==par2)
    autoCorr = true;
else
    autoCorr = false;
end;
autoCorr = false;
% Parameters of the simulation (in the order providen by kernExtractParam)

m = [lfmKern1.mass lfmKern2.mass];                  % Par. 1
D = [lfmKern1.spring lfmKern2.spring];              % Par. 2
C = [lfmKern1.damper lfmKern2.damper];              % Par. 3
sigma2 = 2/lfmKern1.inverseWidth;                   % Par. 4
sigma = sqrt(sigma2);
S = [lfmKern1.sensitivity lfmKern2.sensitivity];    % Par. 5

alpha = C./(2*m);
omega = sqrt(D./m-alpha.^2);

% Initialization of vectors and matrices

g1 = zeros(1,5);
g2 = zeros(1,5);


% Precomputations 

if isreal(omega)
    gamma1 = alpha(1) + j*omega(1);
    gamma2 = alpha(2) + j*omega(2);
    ComputeH1 = lfmComputeH(conj(gamma2),gamma1,sigma2,t1,t2,1);
    ComputeH2 = lfmComputeH(gamma1,conj(gamma2),sigma2,t2,t1,2);
    ComputeH3 = lfmComputeH(gamma2,gamma1,sigma2,t1,t2,1);
    ComputeH4 = lfmComputeH(gamma1,gamma2,sigma2,t2,t1,2);    
    preKernel = ComputeH1 + ComputeH2 - ComputeH3 -ComputeH4;    
    ComputeUpsilon1 = lfmComputeUpsilon(conj(gamma2),sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon2 = lfmComputeUpsilon(gamma1,sigma2,t1,zeros(size(t2)),4);    
    ComputeUpsilon3 = lfmComputeUpsilon(gamma2,sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon4 = lfmComputeUpsilon(gamma1,sigma2,t1,zeros(size(t2)),4);
else
    gamma1_p = alpha(1) + j*omega(1);
    gamma1_m = alpha(1) - j*omega(1);
    gamma2_p = alpha(2) + j*omega(2);
    gamma2_m = alpha(2) - j*omega(2);
    ComputeH1 = lfmComputeH(gamma2_m,gamma1_p,sigma2,t1,t2,1);
    ComputeH2 = lfmComputeH(gamma1_p,gamma2_m,sigma2,t2,t1,2);
    ComputeH3 = lfmComputeH(gamma2_p,gamma1_m,sigma2,t1,t2,1);
    ComputeH4 = lfmComputeH(gamma1_m,gamma2_p,sigma2,t2,t1,2);
    ComputeH5 = lfmComputeH(gamma2_m,gamma1_m,sigma2,t1,t2,1);
    ComputeH6 = lfmComputeH(gamma1_m,gamma2_m,sigma2,t2,t1,2);
    ComputeH7 = lfmComputeH(gamma2_p,gamma1_p,sigma2,t1,t2,1);
    ComputeH8 = lfmComputeH(gamma1_p,gamma2_p,sigma2,t2,t1,2);
    preKernel = ComputeH1 + ComputeH2 + ComputeH3 + ComputeH4 ...
        - ComputeH5 - ComputeH6 - ComputeH7 - ComputeH8;
    ComputeUpsilon1 = lfmComputeUpsilon(gamma2_m,sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon2 = lfmComputeUpsilon(gamma1_p,sigma2,t1,zeros(size(t2)),4);
    ComputeUpsilon3 = lfmComputeUpsilon(gamma2_p,sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon4 = lfmComputeUpsilon(gamma1_m,sigma2,t1,zeros(size(t2)),4);
    ComputeUpsilon5 = lfmComputeUpsilon(gamma2_m,sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon6 = lfmComputeUpsilon(gamma1_m,sigma2,t1,zeros(size(t2)),4);
    ComputeUpsilon7 = lfmComputeUpsilon(gamma2_p,sigma2,t2,zeros(size(t1)),3);
    ComputeUpsilon8 = lfmComputeUpsilon(gamma1_p,sigma2,t1,zeros(size(t2)),4);
end


% Gradient with respect to m, D and C
for ind_theta = 1:3 % Parameter (m, D or C)
  for ind_par = 0:1 % System (1 or 2)

    % Choosing the right gradients for m, omega, gamma1 and gamma2

    switch ind_theta
        case 1  % Gradient wrt m
            gradThetaM = [1-ind_par ind_par];
            gradThetaAlpha = -C./(2*(m.^2));
            gradThetaOmega = (C.^2-2*m.*D)./(2*(m.^2).*sqrt(4*m.*D-C.^2));
        case 2  % Gradient wrt D
            gradThetaM = zeros(1,2);
            gradThetaAlpha = zeros(1,2);
            gradThetaOmega = 1./sqrt(4*m.*D-C.^2);
        case 3  % Gradient wrt C
            gradThetaM = zeros(1,2);
            gradThetaAlpha = 1./(2*m);
            gradThetaOmega = -C./(2*m.*sqrt(4*m.*D-C.^2));
    end
  
    gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
    gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
  
    % Gradient evaluation
    
    if isreal(omega)
        gamma1 = alpha(1) + j*omega(1);
        gamma2 = alpha(2) + j*omega(2);

        if autoCorr
            matInd = ones(4,2);
        else
            matInd = repmat([ind_par 1-ind_par; 1-ind_par ind_par],2,1);
        end;
        
        gradThetaGamma = [conj(gradThetaGamma1(2)), gradThetaGamma1(1); ...
            gradThetaGamma1(1), conj(gradThetaGamma1(2)); ...
            gradThetaGamma1(2), gradThetaGamma1(1); ...
            gradThetaGamma1(1), gradThetaGamma1(2)].*matInd;
        matGrad = (sigma*prod(S)*sqrt(pi)/(4*prod(m)*prod(omega))) ...
                * real(lfmGradientH(conj(gamma2),gamma1,sigma2, ...
                    gradThetaGamma(1,:),t1, t2, 1, ComputeH1, ComputeUpsilon1) ...
                + lfmGradientH(gamma1,conj(gamma2),sigma2, ...
                    gradThetaGamma(2,:),t2, t1, 2, ComputeH2, ComputeUpsilon2) ...
                - lfmGradientH(gamma2,gamma1,sigma2, ...
                    gradThetaGamma(3,:),t1, t2, 1, ComputeH3, ComputeUpsilon3) ...
                - lfmGradientH(gamma1,gamma2,sigma2, ...
                    gradThetaGamma(4,:),t2, t1, 2, ComputeH4, ComputeUpsilon4) ...
                - (1+autoCorr)*(gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);
    else
        gamma1_p = alpha(1) + j*omega(1);
        gamma1_m = alpha(1) - j*omega(1);
        gamma2_p = alpha(2) + j*omega(2);
        gamma2_m = alpha(2) - j*omega(2);
        
        if autoCorr
            matInd = ones(8,2);
        else
            matInd = repmat([ind_par 1-ind_par; 1-ind_par ind_par],4,1);
        end;

        gradThetaGamma = [gradThetaGamma2(2) gradThetaGamma1(1); ...
            gradThetaGamma1(1) gradThetaGamma2(2); ...
            gradThetaGamma1(2) gradThetaGamma2(1); ...
            gradThetaGamma2(1) gradThetaGamma1(2); ...
            gradThetaGamma2(2) gradThetaGamma2(1); ...
            gradThetaGamma2(1) gradThetaGamma2(2); ...
            gradThetaGamma1(2) gradThetaGamma1(1); ...
            gradThetaGamma1(1) gradThetaGamma1(2)].*matInd;
        matGrad = (sigma*prod(S)*sqrt(pi)/(8*prod(m)*prod(omega))) ...
                * (lfmGradientH(gamma2_m,gamma1_p,sigma2, ...
                    gradThetaGamma(1,:),t1, t2, 1, ComputeH1, ComputeUpsilon1) ...
                + lfmGradientH(gamma1_p,gamma2_m,sigma2, ...
                    gradThetaGamma(2,:),t2,t1, 2, ComputeH2, ComputeUpsilon2) ...
                + lfmGradientH(gamma2_p,gamma1_m,sigma2, ...
                    gradThetaGamma(3,:),t1, t2, 1, ComputeH3, ComputeUpsilon3) ...
                + lfmGradientH(gamma1_m,gamma2_p,sigma2, ...
                    gradThetaGamma(4,:),t2 ,t1 , 2, ComputeH4, ComputeUpsilon4) ...
                - lfmGradientH(gamma2_m,gamma1_m,sigma2, ...
                    gradThetaGamma(5,:),t1, t2, 1, ComputeH5, ComputeUpsilon5) ...
                - lfmGradientH(gamma1_m,gamma2_m,sigma2, ...
                    gradThetaGamma(6,:),t2, t1, 2, ComputeH6, ComputeUpsilon6) ...
                - lfmGradientH(gamma2_p,gamma1_p,sigma2, ...
                    gradThetaGamma(7,:),t1, t2, 1, ComputeH7, ComputeUpsilon7) ...
                - lfmGradientH(gamma1_p,gamma2_p,sigma2, ...
                    gradThetaGamma(8,:),t2, t1, 2, ComputeH8, ComputeUpsilon8) ...
                - (1+autoCorr)*(gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);
    end
    % Check the parameter to assign the derivative
    if ind_par == 0
      g1(ind_theta) = sum(sum(matGrad.*covGrad));
    else
      g2(ind_theta) = sum(sum(matGrad.*covGrad));
    end
    % If we are evaluating an autocorrelation there is no need to obtain
    % the second derivative
    if autoCorr
        g2(ind_theta) = g1(ind_theta);
        break;
    end;
  end
end

% Gradients with respect to sigma

if isreal(omega)
    gamma1 = alpha(1) + j*omega(1);
    gamma2 = alpha(2) + j*omega(2);
    matGrad = (prod(S)*sqrt(pi)/(4*prod(m)*prod(omega))) ...
        * real(preKernel ...
        + sigma*(lfmGradientSigmaH(conj(gamma2),gamma1,sigma2,t1,t2,1) ...
        + lfmGradientSigmaH(gamma1,conj(gamma2),sigma2,t2,t1,2) ...
        - lfmGradientSigmaH(gamma2,gamma1,sigma2,t1,t2,1) ...
        - lfmGradientSigmaH(gamma1,gamma2,sigma2,t2,t1,2)));
else
    gamma1_p = alpha(1) + j*omega(1);
    gamma1_m = alpha(1) - j*omega(1);
    gamma2_p = alpha(2) + j*omega(2);
    gamma2_m = alpha(2) - j*omega(2);
    matGrad = (prod(S)*sqrt(pi)/(8*prod(m)*prod(omega))) ...
       * (preKernel ...
       + sigma*(lfmGradientSigmaH(gamma2_m,gamma1_p,sigma2,t1,t2, 1) ...
       + lfmGradientSigmaH(gamma1_p,gamma2_m,sigma2,t2,t1, 2) ...
       + lfmGradientSigmaH(gamma2_p,gamma1_m,sigma2,t1,t2, 1) ...
       + lfmGradientSigmaH(gamma1_m,gamma2_p,sigma2,t2,t1, 2) ...
       - lfmGradientSigmaH(gamma2_m,gamma1_m,sigma2,t1,t2, 1) ...
       - lfmGradientSigmaH(gamma1_m,gamma2_m,sigma2,t2,t1, 2) ...
       - lfmGradientSigmaH(gamma2_p,gamma1_p,sigma2,t1,t2, 1) ...
       - lfmGradientSigmaH(gamma1_p,gamma2_p,sigma2,t2,t1, 2)));
end;

g1(4) = sum(sum(matGrad.*covGrad))*(-(sigma^3)/4);
g2(4) = g1(4);

% Gradients with respect to S

if isreal(omega)
    matGrad = (sigma*sqrt(pi)/(4*prod(m)*prod(omega))) ...
            * real(preKernel);
else
    matGrad = (sigma*sqrt(pi)/(8*prod(m)*prod(omega))) ...
           * (preKernel);
end;

g1(5) = sum(sum((1+autoCorr)*S(2)*matGrad.*covGrad));
if autoCorr
    g2(5) = g1(5);
else
    g2(5) = sum(sum(S(1)*matGrad.*covGrad));
end
g2(4) = 0; % Otherwise is counted twice, temporarly changed by Mauricio Alvarez 
return