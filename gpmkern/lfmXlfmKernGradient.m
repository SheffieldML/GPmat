function [g1, g2] = lfmXlfmKernGradient(lfmKern1, lfmKern2, t1, t2, covGrad, meanVector)

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
% ARG meanVec : precomputed factor that is used for the switching dynamical
% latent force model.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see lfmKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see lfmKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, lfmKernParamInit, lfmKernExtractParam
%
% COPYRIGHT : David Luengo, 2007, 2008, Mauricio Alvarez, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010

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

if isreal(omega)
    computeH = cell(4,1);
    computeUpsilonMatrix = cell(2,1);
    computeUpsilonVector = cell(2,1);
    gradientUpsilonMatrix = cell(2,1);
    gradientUpsilonVector = cell(2,1);
    gamma1 = alpha(1) + j*omega(1);
    gamma2 = alpha(2) + j*omega(2);
    gradientUpsilonMatrix{1} = lfmGradientUpsilonMatrix(gamma1,sigma2,t1, t2);
    gradientUpsilonMatrix{2} = lfmGradientUpsilonMatrix(gamma2,sigma2,t2, t1);
    gradientUpsilonVector{1} = lfmGradientUpsilonVector(gamma1,sigma2,t1);
    gradientUpsilonVector{2} = lfmGradientUpsilonVector(gamma2,sigma2,t2);
    preGamma(1) = gamma1 + gamma2;
    preGamma(2) = conj(gamma1) + gamma2;
    preGamma2 = preGamma.^2;
    preConst = 1./preGamma;
    preConst2 = 1./preGamma2;
    preExp1 = exp(-gamma1*t1);
    preExp2 = exp(-gamma2*t2);
    preExpt1 = t1.*exp(-gamma1*t1);
    preExpt2 = t2.*exp(-gamma2*t2);
    [computeH{1}, computeUpsilonMatrix{1}] = lfmComputeH3(gamma1, gamma2, sigma2, t1,t2,preConst, 0, 1);
    [computeH{2}, computeUpsilonMatrix{2}] = lfmComputeH3(gamma2, gamma1, sigma2, t2,t1,preConst(2) - preConst(1), 0, 0);
    [computeH{3}, computeUpsilonVector{1}] = lfmComputeH4(gamma1, gamma2, sigma2, t1, preGamma, preExp2, 0, 1  );
    [computeH{4}, computeUpsilonVector{2}] = lfmComputeH4(gamma2, gamma1, sigma2, t2, preGamma, preExp1,0, 0 );
    preKernel = real( computeH{1} + computeH{2}.' + computeH{3} + computeH{4}.');
else
    computeH  = cell(4,1);
    computeUpsilonMatrix = cell(2,1);
    computeUpsilonVector = cell(2,1);
    gradientUpsilonMatrix = cell(4,1);
    gradientUpsilonVector = cell(4,1);
    gamma1_p = alpha(1) + j*omega(1);
    gamma1_m = alpha(1) - j*omega(1);
    gamma2_p = alpha(2) + j*omega(2);
    gamma2_m = alpha(2) - j*omega(2);
    gradientUpsilonMatrix{1} = lfmGradientUpsilonMatrix(gamma1_p,sigma2,t1, t2);
    gradientUpsilonMatrix{2} = lfmGradientUpsilonMatrix(gamma1_m,sigma2,t1, t2);
    gradientUpsilonMatrix{3} = lfmGradientUpsilonMatrix(gamma2_p,sigma2,t2, t1);
    gradientUpsilonMatrix{4} = lfmGradientUpsilonMatrix(gamma2_m,sigma2,t2, t1);
    gradientUpsilonVector{1} = lfmGradientUpsilonVector(gamma1_p,sigma2,t1);
    gradientUpsilonVector{2} = lfmGradientUpsilonVector(gamma1_m,sigma2,t1);
    gradientUpsilonVector{3} = lfmGradientUpsilonVector(gamma2_p,sigma2,t2);
    gradientUpsilonVector{4} = lfmGradientUpsilonVector(gamma2_m,sigma2,t2);
    preExp1 = zeros(length(t1),2);
    preExp2 = zeros(length(t2),2);
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
    preExp2(:,1) = exp(-gamma2_p*t2);
    preExp2(:,2) = exp(-gamma2_m*t2);
    preExpt1(:,1) = t1.*exp(-gamma1_p*t1);
    preExpt1(:,2) = t1.*exp(-gamma1_m*t1);
    preExpt2(:,1) = t2.*exp(-gamma2_p*t2);
    preExpt2(:,2) = t2.*exp(-gamma2_m*t2);
    [computeH{1}, computeUpsilonMatrix{1}] =  lfmComputeH3(gamma1_p, gamma1_m, sigma2, t1,t2,preFactors([1 2]), 1);
    [computeH{2}, computeUpsilonMatrix{2}] =  lfmComputeH3(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1);
    [computeH{3}, computeUpsilonVector{1}] =  lfmComputeH4(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 1 );
    [computeH{4}, computeUpsilonVector{2}] =  lfmComputeH4(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 );
    preKernel = ( computeH{1} + computeH{2}.' + computeH{3} + computeH{4}.');
end

if isreal(omega)
    if lfmKern1.isNormalised
        K0 = (prod(S)/(4*sqrt(2)*prod(m)*prod(omega)));
    else
        K0 = (sigma*prod(S)*sqrt(pi)/(4*prod(m)*prod(omega)));
    end
else    
    if lfmKern1.isNormalised
        K0 = (prod(S)/(8*sqrt(2)*prod(m)*prod(omega)));
    else
        K0 = (sigma*prod(S)*sqrt(pi)/(8*prod(m)*prod(omega)));
    end
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
            
            gradThetaGamma2 = gradThetaGamma1(2);
            gradThetaGamma1 = [gradThetaGamma1(1) conj(gradThetaGamma1(1))];
            
            %  gradThetaGamma1 =  gradThetaGamma11;
            
            
            if ~ind_par
                matGrad = K0* ...
                     real( lfmGradientH31( preConst, preConst2, gradThetaGamma1, ...
                    gradientUpsilonMatrix{1}, 1, computeUpsilonMatrix{1}, 1, 0, 1) + ...
                    lfmGradientH32( preGamma2, gradThetaGamma1, computeUpsilonMatrix{2}, ...
                    1, 0, 0).' + ...
                    lfmGradientH41( preGamma, preGamma2, gradThetaGamma1, preExp2, ...
                    gradientUpsilonVector{1}, 1, computeUpsilonVector{1}, 1, 0, 1) + ...
                    lfmGradientH42(preGamma, preGamma2, gradThetaGamma1, preExp1, preExpt1, ...
                    computeUpsilonVector{2}, 1, 0, 0).' ...
                    - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                    + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                    *preKernel);
            else
                matGrad = K0* ...
                     real( lfmGradientH31( (preConst(2)- preConst(1)), (-preConst2(2) + preConst2(1)), gradThetaGamma2, ...
                    gradientUpsilonMatrix{2}, 1, computeUpsilonMatrix{2}, 1, 0, 0).' + ...
                    lfmGradientH32( preConst2, gradThetaGamma2, computeUpsilonMatrix{1}, ...
                    1, 0, 1) + ...
                    lfmGradientH41( preGamma, preGamma2, gradThetaGamma2, preExp1, ...
                    gradientUpsilonVector{2}, 1, computeUpsilonVector{2}, 1, 0, 0).' + ...
                    lfmGradientH42(preGamma, preGamma2, gradThetaGamma2, preExp2, preExpt2, ...
                    computeUpsilonVector{1}, 1, 0, 1) ...
                    - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                    + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                    *preKernel);
            end
        else
                     
            gradThetaGamma11 = [gradThetaGamma1(1) gradThetaGamma2(1)];
            gradThetaGamma2 = [gradThetaGamma1(2) gradThetaGamma2(2)];
            gradThetaGamma1 = gradThetaGamma11;
            
            if ~ind_par % ind_par = k
                matGrad = K0* ...
                     ( lfmGradientH31( preFactors([1 2]), preFactors2([1 2]), gradThetaGamma1, ...
                    gradientUpsilonMatrix{1}, gradientUpsilonMatrix{2}, computeUpsilonMatrix{1}{1}, ...
                    computeUpsilonMatrix{1}{2}, 1) + ...
                    lfmGradientH32( preGamma2,  gradThetaGamma1, computeUpsilonMatrix{2}{1}, ...
                    computeUpsilonMatrix{2}{2}, 1).' + ...
                    lfmGradientH41( preGamma, preGamma2, gradThetaGamma1, preExp2, ...
                    gradientUpsilonVector{1}, gradientUpsilonVector{2}, computeUpsilonVector{1}{1},...
                    computeUpsilonVector{1}{2}, 1) + ...
                    lfmGradientH42(preGamma, preGamma2, gradThetaGamma1, preExp1, preExpt1, ...
                    computeUpsilonVector{2}{1}, computeUpsilonVector{2}{2}, 1).'...
                    - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                    + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                    *preKernel);
                
            else % ind_par = r
                matGrad = K0* ...
                     ( lfmGradientH31( preFactors([3 4]), preFactors2([3 4]), gradThetaGamma2, ...
                    gradientUpsilonMatrix{3}, gradientUpsilonMatrix{4}, computeUpsilonMatrix{2}{1}, ...
                    computeUpsilonMatrix{2}{2}, 1).' + ...
                    lfmGradientH32( preGamma2([1 3 2 4]),  gradThetaGamma2, computeUpsilonMatrix{1}{1}, ...
                    computeUpsilonMatrix{1}{2}, 1) + ...
                    lfmGradientH41( preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExp1, ...
                    gradientUpsilonVector{3}, gradientUpsilonVector{4}, computeUpsilonVector{2}{1},...
                    computeUpsilonVector{2}{2}, 1).' + ...
                    lfmGradientH42(preGamma([1 3 2 4]), preGamma2([1 3 2 4]), gradThetaGamma2, preExp2, preExpt2, ...
                    computeUpsilonVector{1}{1}, computeUpsilonVector{1}{2}, 1)...
                    - (gradThetaM(1+ind_par)/m(1+ind_par) ...
                    + gradThetaOmega(1+ind_par)/omega(1+ind_par)) ...
                    *preKernel);
            end
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

if isreal(omega)
    if lfmKern1.isNormalised
        matGrad = K0*real(lfmGradientSigmaH3(gamma1, gamma2, sigma2, t1,t2,preConst, 0, 1)  ...
            +  lfmGradientSigmaH3(gamma2, gamma1, sigma2, t2,t1,preConst(2) - preConst(1), 0, 0).'...
            +  lfmGradientSigmaH4(gamma1, gamma2, sigma2, t1, preGamma, preExp2, 0, 1  )...
            +  lfmGradientSigmaH4(gamma2, gamma1, sigma2, t2, preGamma, preExp1, 0, 0 ).');        
    else
        matGrad = (prod(S)*sqrt(pi)/(4*prod(m)*prod(omega))) ...
            * real(preKernel ...
            + sigma*(lfmGradientSigmaH3(gamma1, gamma2, sigma2, t1,t2,preConst, 0, 1)  ...
            +  lfmGradientSigmaH3(gamma2, gamma1, sigma2, t2,t1,preConst(2) - preConst(1), 0, 0).'...
            +  lfmGradientSigmaH4(gamma1, gamma2, sigma2, t1, preGamma, preExp2, 0, 1  )...
            +  lfmGradientSigmaH4(gamma2, gamma1, sigma2, t2, preGamma, preExp1, 0, 0 ).'));
    end
else
    if lfmKern1.isNormalised
        matGrad = K0*(lfmGradientSigmaH3(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 1)...
            +  lfmGradientSigmaH3(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
            +  lfmGradientSigmaH4(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 1 )...
            +  lfmGradientSigmaH4(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 ).' );       
    else
        matGrad = (prod(S)*sqrt(pi)/(8*prod(m)*prod(omega))) ...
            * (preKernel ...
            + sigma*(lfmGradientSigmaH3(gamma1_p, gamma1_m, sigma2, t1, t2, preFactors([1 2]), 1)...
            +  lfmGradientSigmaH3(gamma2_p, gamma2_m, sigma2, t2, t1, preFactors([3 4]), 1).'...
            +  lfmGradientSigmaH4(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 1 )...
            +  lfmGradientSigmaH4(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 ).' ));
    end
end;

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

if isreal(omega)
    if lfmKern1.isNormalised
        matGrad = (1/(4*sqrt(2)*prod(m)*prod(omega))) ...
            * real(preKernel);
    else
        matGrad = (sigma*sqrt(pi)/(4*prod(m)*prod(omega))) ...
            * real(preKernel);
    end
else
    if lfmKern1.isNormalised
        matGrad = (1/(8*sqrt(2)*prod(m)*prod(omega))) ...
            * (preKernel);
    else
        matGrad = (sigma*sqrt(pi)/(8*prod(m)*prod(omega))) ...
            * (preKernel);
    end
end;

if subComponent
    if size(meanVector,1) ==1,
        matGrad = matGrad*meanVector;
    else
        matGrad = (meanVector*matGrad).';
    end
end

g1(5) = sum(sum(S(2)*matGrad.*covGrad));
g2(5) = sum(sum(S(1)*matGrad.*covGrad));

g2(4) = 0; % Otherwise is counted twice

g1 = real(g1);
g2 = real(g2);

return
