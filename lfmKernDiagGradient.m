function g = lfmKernDiagGradient(lfmKern, t, covDiag)

% LFMKERNDIAGGRADIENT Compute the gradient of the LFM kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% single input motif kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% lfmKernExtractParam command.
% ARG lfmKern : the kernel structure for which the gradients are
% computed.
% ARG t : the input data for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% lfmKernExtractParam.
%
% SEEALSO : lfmKernParamInit, kernDiagGradient, lfmKernExtractParam, lfmKernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(t, 2) > 1
    error('Input can only have one column');
end

% Parameters of the simulation (in the order providen by kernExtractParam)
m = lfmKern.mass;                  % Par. 1
D = lfmKern.spring;                % Par. 2
C = lfmKern.damper;                % Par. 3
sigma2 = 2/lfmKern.inverseWidth;   % Par. 4
sigma = sqrt(sigma2);
S = lfmKern.sensitivity ;          % Par. 5

alpha = C/(2*m);
omega = sqrt(D/m-alpha^2);

% Initialization of vectors and matrices
g = zeros(1,5);

% Precomputations
diagH = cell(1,4);
gradDiag = cell(1,2);
preExp = zeros(length(t),2);
gamma_p = alpha + j*omega;
gamma_m = alpha - j*omega;
preFactors(1) = 2/(gamma_p + gamma_m) - 1/gamma_m;
preFactors(2) = 2/(gamma_p + gamma_m) - 1/gamma_p;

preExp(:,1) = exp(-gamma_p*t);
preExp(:,2) = exp(-gamma_m*t);
% Actual computation of the kernel
[diagH{1}, upsilonDiag{2}]  = lfmDiagComputeH3(-gamma_m, sigma2, t, preFactors(1), preExp(:,2), 1);
[diagH{2}, upsilonDiag{1}]  = lfmDiagComputeH3(-gamma_p, sigma2, t, preFactors(2), preExp(:,1), 1);
[diagH{3}, upsilonDiag{4}]  = lfmDiagComputeH4( gamma_m, sigma2, t, [gamma_m  (gamma_p + gamma_m)], [preExp(:,2) preExp(:,1)] , 1);
[diagH{4}, upsilonDiag{3}]  = lfmDiagComputeH4( gamma_p, sigma2, t, [gamma_p  (gamma_p + gamma_m)], preExp , 1);


gradDiag{1} =  lfmGradientUpsilonVector(-gamma_p, sigma2, t);
gradDiag{2} =  lfmGradientUpsilonVector(-gamma_m, sigma2, t);
gradDiag{3} =  lfmGradientUpsilonVector(gamma_p, sigma2, t);
gradDiag{4} =  lfmGradientUpsilonVector(gamma_m, sigma2, t);

preKernel = diagH{1} + diagH{2} + diagH{3} + diagH{4};

if lfmKern.isNormalised
    k0 = lfmKern.sensitivity^2/(8*sqrt(2)*lfmKern.mass^2*omega^2);
else
    k0 = sqrt(pi)*sigma*lfmKern.sensitivity^2/(8*lfmKern.mass^2*omega^2);
end


% Gradient with respect to m, D and C
for ind_theta = 1:3 % Parameter (m, D or C)
    % Choosing the right gradients for m, omega, gamma1 and gamma2
    switch ind_theta
        case 1  % Gradient wrt m
            gradThetaM = 1;
            gradThetaAlpha = -C/(2*(m^2));
            gradThetaOmega = (C^2-2*m*D)/(2*(m^2)*sqrt(4*m*D-C^2));
        case 2  % Gradient wrt D
            gradThetaM = 0;
            gradThetaAlpha = zeros(1,2);
            gradThetaOmega = 1/sqrt(4*m*D-C^2);
        case 3  % Gradient wrt C
            gradThetaM = 0;
            gradThetaAlpha = 1/(2*m);
            gradThetaOmega = -C/(2*m*sqrt(4*m*D-C^2));
    end

    gradThetaGamma1 = gradThetaAlpha + j*gradThetaOmega;
    gradThetaGamma2 = gradThetaAlpha - j*gradThetaOmega;
    % Gradient evaluation
    gradThetaGamma = [gradThetaGamma1(1) gradThetaGamma2(1)];
    matGrad = lfmDiagGradientH3(- gamma_m, t,preFactors(1), preExp(:,2), ...
        upsilonDiag{2}, gradDiag{2}, diagH{1}, gamma_p + gamma_m, gradThetaGamma) ...
        + lfmDiagGradientH3(- gamma_p, t,preFactors(2), preExp(:,1), ...
        upsilonDiag{1}, gradDiag{1}, diagH{2}, gamma_p + gamma_m, [gradThetaGamma(2) gradThetaGamma(1)]) ...
        + lfmDiagGradientH4( t, [gamma_m  (gamma_p + gamma_m)], [preExp(:,2) preExp(:,1)], ...
        upsilonDiag{4}, gradDiag{4}, gradThetaGamma) ...
        + lfmDiagGradientH4( t, [gamma_p  (gamma_p + gamma_m)], preExp, ...
        upsilonDiag{3}, gradDiag{3}, [gradThetaGamma(2) gradThetaGamma(1)]) ...
        - 2*(gradThetaM/m + gradThetaOmega/omega)*preKernel;
    g(ind_theta) = k0*sum(sum(matGrad.*covDiag));
end

% Gradients with respect to sigma
if lfmKern.isNormalised
    matGrad = k0*(lfmDiagGradientSH3(-gamma_m, sigma2, t, preFactors(1), preExp(:,2), 1) + ...
        lfmDiagGradientSH3(-gamma_p, sigma2, t, preFactors(2), preExp(:,1), 1) + ...
        lfmDiagGradientSH4( gamma_m, sigma2, t, [gamma_m  (gamma_p + gamma_m)], [preExp(:,2) preExp(:,1)] , 1) + ...
        lfmDiagGradientSH4( gamma_p, sigma2, t, [gamma_p  (gamma_p + gamma_m)], preExp , 1));
else
    matGrad = (S^2*sqrt(pi)/(8*m^2*omega^2)) ...
        * (preKernel ...
        + sigma*(lfmDiagGradientSH3(-gamma_m, sigma2, t, preFactors(1), preExp(:,2), 1) + ...
        lfmDiagGradientSH3(-gamma_p, sigma2, t, preFactors(2), preExp(:,1), 1) + ...
        lfmDiagGradientSH4( gamma_m, sigma2, t, [gamma_m  (gamma_p + gamma_m)], [preExp(:,2) preExp(:,1)] , 1) + ...
        lfmDiagGradientSH4( gamma_p, sigma2, t, [gamma_p  (gamma_p + gamma_m)], preExp , 1)));
end

g(4) = sum(sum(matGrad.*covDiag))*(-(sigma^3)/4);

% Gradients with respect to S

if lfmKern.isNormalised
    matGrad = (1/(8*sqrt(2)*m^2*omega^2)) ...
        * (preKernel);
else
    matGrad = (sigma*sqrt(pi)/(8*m^2*omega^2)) ...
        * (preKernel);
end
g(5) = 2*S*sum(sum(matGrad.*covDiag));

g = real(g);

end

function [vec, upsi] = lfmDiagComputeH3(gamma, sigma2, t, factor, preExp, mode)

upsi = lfmComputeUpsilonVector(gamma ,sigma2, t);
if mode
    vec = preExp.*upsi*factor;
else
    temp = preExp.*upsi;
    vec = 2*real(temp/factor(1)) - temp/factor(2);
end
end

function [vec, upsi] = lfmDiagComputeH4(gamma, sigma2, t, factor, preExp, mode)
upsi = lfmComputeUpsilonVector(gamma ,sigma2, t);
if mode
    vec = (preExp(:,1)/factor(1) -  2*preExp(:,2)/factor(2)).*upsi;
else
    temp2 = upsi.*conj(preExp)/factor(2);
    vec = upsi.*preExp/factor(1) - 2*real(temp2);
end
end

function pgrad = lfmDiagGradientH3(gamma, t, factor, preExp, ...
    compUpsilon, gradUpsilon, termH, preFactorGrad, gradTheta)

expUpsilon = preExp.*compUpsilon;

pgrad = - expUpsilon*(2/preFactorGrad^2)*gradTheta(1) - (t.*termH ...
    + preExp.*gradUpsilon*factor(1) - expUpsilon*(1/gamma^2 - 2/preFactorGrad^2))*gradTheta(2);
end

function pgrad = lfmDiagGradientH4( t, factor, preExp, ...
    compUpsilon, gradUpsilon, gradTheta)

pgrad = 2*preExp(:,2).*compUpsilon.*(t/factor(2) + 1/factor(2)^2)*gradTheta(1) ...
    + (gradUpsilon.*(preExp(:,1)/factor(1) - 2*preExp(:,2)/factor(2)) ...
    - compUpsilon.*(t.*preExp(:,1)/factor(1) + preExp(:,1)/factor(1)^2 ...
    -  2*preExp(:,2)/factor(2)^2))*gradTheta(2); 

end

function [vec, upsi] = lfmDiagGradientSH3(gamma, sigma2, t, factor, preExp, mode)

upsi = lfmGradientSigmaUpsilonVector(gamma ,sigma2, t);
if mode
    vec = preExp.*upsi*factor;
else
    temp = preExp.*upsi;
    vec = 2*real(temp/factor(1)) - temp/factor(2);
end
end

function [vec, upsi] = lfmDiagGradientSH4(gamma, sigma2, t, factor, preExp, mode)

upsi = lfmGradientSigmaUpsilonVector(gamma ,sigma2, t);

if mode
    vec = (preExp(:,1)/factor(1) -  2*preExp(:,2)/factor(2)).*upsi;
else
    temp2 = upsi.*conj(preExp)/factor(2);
    vec = upsi.*preExp/factor(1) - 2*real(temp2);
end
end
