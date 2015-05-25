function [g1, g2] = lfmvXlfmKernGradient(lfmKern1, lfmKern2, t1, t2, covGrad, meanVector)

% LFMVXLFMKERNGRADIENT Compute a cross gradient for a LFMVXLFM.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two lfm kernels, one lfm kernel corresponds to the velocity and
% the other corresponds to the position. It is supposed to be used together
% with the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel (velocity).
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel (position).
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two lfm kernels, one lfm kernel corresponds to the velocity and
% the other corresponds to the position. It is supposed to be used together
% with the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel (velocity).
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel (position).
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two lfm kernels, one lfm kernel corresponds to the velocity and
% the other corresponds to the position. It is supposed to be used together
% with the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel (velocity).
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel (position).
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG meanVec : precomputed factor that is used for the switching dynamical
% latent force model.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmvKernParamInit
%
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

subComponent = false; % This is just a flag that indicates if this kernel is part of a bigger kernel (SDLFM)

if nargin == 4
    covGrad = t2;
    t2 = t1;
elseif nargin == 6
    subComponent = true;
    if numel(meanVector)>1
        if size(meanVector,1) == 1,
            if size(meanVector, 2)~=size(covGrad, 2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        else
            if size((meanVector'), 2)~=size(covGrad,2)
                error('The dimensions of meanVector don''t correspond to the dimensions of covGrad')
            end
        end
    else
        if numel(t1)==1 && numel(t2)>1
            % matGrad will be row vector and so should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [1 dimcovGrad]);
        elseif numel(t1)>1 && numel(t2)==1
            % matGrad will be column vector and sp should be covGrad
            dimcovGrad = length(covGrad);
            covGrad = reshape(covGrad, [dimcovGrad 1]);
        end
    end
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

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
computeH  = cell(4,1);
computeUpsilonMatrix = cell(2,1);
computeUpsilonVector = cell(2,1);
gradientUpsilonMatrix = cell(4,1);
gradientUpsilonVector = cell(4,1);
gamma1_p = alpha(1) + j*omega(1);
gamma1_m = alpha(1) - j*omega(1);
gamma2_p = alpha(2) + j*omega(2);
gamma2_m = alpha(2) - j*omega(2);
preExp1 = zeros(length(t1),2);
preExp2 = zeros(length(t2),2);
preExpg1 = zeros(length(t1),2);
preExpt1 = zeros(length(t1),2);
preExpt2 = zeros(length(t2),2);
preGamma(1) = gamma1_p + gamma2_p;
preGamma(2) = gamma1_p + gamma2_m;
preGamma(3) = gamma1_m + gamma2_p;
preGamma(4) = gamma1_m + gamma2_m;
preGamma2 = preGamma.^2;
preConst = 1./preGamma;
preConst2 = 1./(preGamma2);
preFactors(1) = preConst(2) - preConst(1);
preFactors(2) = preConst(3) - preConst(4);
preFactors(3) = preConst(3) - preConst(1);
preFactors(4) = preConst(2) - preConst(4);
preFactors2(1) = -preConst2(2) + preConst2(1);
preFactors2(2) = -preConst2(3) + preConst2(4);
preFactors2(3) = -preConst2(3) + preConst2(1);
preFactors2(4) = -preConst2(2) + preConst2(4);
preExp1(:,1) = exp(-gamma1_p*t1);
preExp1(:,2) = exp(-gamma1_m*t1);
preExpg1(:,1) = gamma1_p*preExp1(:,1);
preExpg1(:,2) = gamma1_m*preExp1(:,2);
preExp2(:,1) = exp(-gamma2_p*t2);
preExp2(:,2) = exp(-gamma2_m*t2);
preExpt1(:,1) = t1.*preExpg1(:,1);
preExpt1(:,2) = t1.*preExpg1(:,2);
preExpt2(:,1) = t2.*exp(-gamma2_p*t2);
preExpt2(:,2) = t2.*exp(-gamma2_m*t2);
[computeH{1}, computeUpsilonMatrix{1}, computeUpsilonMatrixLocal{1}] =  lfmComputeH3VP(gamma1_p, gamma1_m, sigma2, t1,t2, preFactors([1 2]), 0);
[computeH{2}, computeUpsilonMatrix{2}, computeUpsilonMatrixLocal{2}] =  lfmComputeH3VP(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1);
[computeH{3}, computeUpsilonVector{1}, computeUpsilonVectorLocal{1}] =  lfmComputeH4VP(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 0 );
[computeH{4}, computeUpsilonVector{2}, computeUpsilonVectorLocal{2}] =  lfmComputeH4VP(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpg1, 1 );
preKernel = ( computeH{1} + computeH{2}.' + computeH{3} + computeH{4}.');

% Precompute derivatives
gradientUpsilonMatrix{1} = lfmvpGradientUpsilonMatrix(gamma1_p,sigma2,t1, t2, 0, computeUpsilonMatrixLocal{1}{1});
gradientUpsilonMatrix{2} = lfmvpGradientUpsilonMatrix(gamma1_m,sigma2,t1, t2, 0, computeUpsilonMatrixLocal{1}{2});
gradientUpsilonMatrix{3} = lfmvpGradientUpsilonMatrix(gamma2_p,sigma2,t2, t1, 1, computeUpsilonMatrixLocal{2}{1});
gradientUpsilonMatrix{4} = lfmvpGradientUpsilonMatrix(gamma2_m,sigma2,t2, t1, 1, computeUpsilonMatrixLocal{2}{2});
gradientUpsilonVector{1} = lfmvpGradientUpsilonVector(gamma1_p,sigma2,t1, computeUpsilonVectorLocal{1}{1});
gradientUpsilonVector{2} = lfmvpGradientUpsilonVector(gamma1_m,sigma2,t1, computeUpsilonVectorLocal{1}{2});
gradientUpsilonVector{3} = lfmGradientUpsilonVector(gamma2_p,sigma2,t2);
gradientUpsilonVector{4} = lfmGradientUpsilonVector(gamma2_m,sigma2,t2);

if lfmKern1.isNormalised
    K0 =  lfmKern1.sensitivity*lfmKern2.sensitivity/(8*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
    K02 = 1/(8*sqrt(2)*prod(m)*prod(omega));
else
    K0 =  sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity/(8*lfmKern1.mass*lfmKern2.mass*prod(omega)); 
    K02 = sigma*sqrt(pi)/(8*prod(m)*prod(omega));
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
        gradThetaGamma11 = [gradThetaGamma1(1) gradThetaGamma2(1)];
        gradThetaGamma2 = [gradThetaGamma1(2) gradThetaGamma2(2)];
        gradThetaGamma1 = gradThetaGamma11;

        if ~ind_par % ind_par = k or d
            matGrad = K0* ( lfmGradientH31( preFactors([1 2]), preFactors2([1 2]), gradThetaGamma1, ...
                gradientUpsilonMatrix{1}, gradientUpsilonMatrix{2}, computeUpsilonMatrix{1}{1}, ...
                computeUpsilonMatrix{1}{2}, 1) + ...
                lfmGradientH32( preGamma2,  gradThetaGamma1, computeUpsilonMatrix{2}{1}, ...
                computeUpsilonMatrix{2}{2}, 1).' + ...
                lfmGradientH41( preGamma, preGamma2, gradThetaGamma1, preExp2, ...
                gradientUpsilonVector{1}, gradientUpsilonVector{2}, computeUpsilonVector{1}{1},...
                computeUpsilonVector{1}{2}, 1) + ...
                lfmGradientH42VP(preGamma, preGamma2, gradThetaGamma1, preExp1, preExpg1, preExpt1, ...
                computeUpsilonVector{2}{1}, computeUpsilonVector{2}{2}).'...
                - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);
        else % ind_par = r or d'
            matGrad = K0*( lfmGradientH31( preFactors([3 4]), preFactors2([3 4]), gradThetaGamma2, ...
                gradientUpsilonMatrix{3}, gradientUpsilonMatrix{4}, computeUpsilonMatrix{2}{1}, ...
                computeUpsilonMatrix{2}{2}, 1).' + ...
                lfmGradientH32( preGamma2([1 3 2 4]),  gradThetaGamma2, computeUpsilonMatrix{1}{1}, ...
                computeUpsilonMatrix{1}{2}, 1) + ...
                lfmGradientH41VP( preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExpg1, ...
                gradientUpsilonVector{3}, gradientUpsilonVector{4}, computeUpsilonVector{2}{1},...
                computeUpsilonVector{2}{2}).' + ...
                lfmGradientH42(preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExp2, preExpt2, ...
                computeUpsilonVector{1}{1}, computeUpsilonVector{1}{2}, 1)...
                - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                *preKernel);            
        end
        if subComponent
            if size(meanVector,1) ==1,
                matGrad = matGrad*meanVector;
            else
                matGrad = (meanVector*matGrad).';
            end
        end
        % Check the parameter to assign the derivative
        if ind_par == 0
            g1(ind_theta) = sum(sum(matGrad.*covGrad));
        else
            g2(ind_theta) = sum(sum(matGrad.*covGrad));
        end       
    end
end

% Gradients with respect to sigma
if lfmKern1.isNormalised    
    matGrad = K0*(lfmGradientSigmaH3VP(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 0)...
        +  lfmGradientSigmaH3VP(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
        +  lfmGradientSigmaH4VP(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 0 )...
        +  lfmGradientSigmaH4VP(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpg1, 1 ).' );    
else
    matGrad = (prod(S)*sqrt(pi)/(8*prod(m)*prod(omega))) ...
        * (preKernel ...
        + sigma*(lfmGradientSigmaH3VP(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 0)...
        +  lfmGradientSigmaH3VP(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
        +  lfmGradientSigmaH4VP(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 0 )...
        +  lfmGradientSigmaH4VP(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExpg1, 1 ).' ));    
end

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(4) = sum(sum(matGrad.*covGrad))*(-(sigma^3)/4);
g2(4) = g1(4);

% Gradients with respect to S
matGrad = K02*preKernel;

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(5) = sum(sum(S(2)*matGrad.*covGrad));
g2(5) = sum(sum(S(1)*matGrad.*covGrad));


g2(4) = 0; % Otherwise is counted twice, temporarly changed by Mauricio Alvarez

g1 = real(g1);
g2 = real(g2);
