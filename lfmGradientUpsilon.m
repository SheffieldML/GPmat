function    g = lfmGradientUpsilon(gamma,sigma2,gradThetaGamma,t1, t2, mode)

% LFMGRADIENTUPSILON Gradient of the function \upsilon(z) with respect to
% one of the hyperparameters of the kernel: m_k, C_k, D_k, m_r, C_r or D_r.
% FORMAT
% DESC Computes the gradient of the function \upsilon(z) with respect to
% one of the parameters of the system (mass, spring or damper).
% ARG gamma : Gamma value of the system.
% ARG sigma2 : length scale of latent process.
% ARG gradThetaGamma : Vector with the gradient of gamma1 and gamma2 with
% respect to the desired parameter.
% ARG Tt1 : first time input (number of time points 1 x number of time points 2).
% ARG Tt2 : second time input (number of time points 1 x number of time points 2).
% ARG mode : indicates in which way the variables t1 and t2 must be
% transposed
% RETURN g : Gradient of the kernel with respect to the desired
% parameter.
%
% COPYRIGHT : David Luengo, Mauricio Alvarez 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientH

% KERN


%Parameters of the function

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

switch mode
    case 1
        % In this mode the Tt1 is actually the Tt1 and Tt2 is actually Tt2
        % Z1 is a matrix with order given by Tt1 and z2 is a vector with
        % order given by Tt2
        Tt11 = repmat(t1, 1, size(t2, 1));
        Tt22 = repmat(t2', size(t1, 1), 1);
        Z1 = (Tt11-Tt22)/sigma - sigma*gamma/2;
        z2 = t2/sigma + sigma*gamma/2;
        Z2 = repmat(z2.', size(t1, 1), 1);
        [preFactorZ1, gradPreFactorZ1] = quadrantWofzGrad(Z1, sigma,1);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma,1);
        Tt1minusTt2 = Tt11 - Tt22;
        preFactor1 = exp(sigma2*(gamma^2)/4 - gamma*(Tt1minusTt2)...
            + log(sigma2*gamma - 2*(Tt1minusTt2)) + log(gradThetaGamma));
        preFactor2 = exp(-(Tt1minusTt2).*(Tt1minusTt2)/sigma2 ...
            + log(gradThetaGamma) + log(gradPreFactorZ1));
        preFactor3 = exp(repmat((-t2.*t2/sigma2).',size(t1, 1), 1) ...
            - gamma*Tt11 + log(gradThetaGamma)...
            + log(Tt11.*repmat(preFactorz2.',size(t1, 1), 1) ...
            + repmat((gradPreFactorz2).',size(t1, 1), 1)));
    case 2

        % In this mode the Tt1 is actually the Tt2 and Tt2 is actually Tt1
        % Z1 is a matrix with the order given by Tt2 and z2 is a vector
        % with the order given by Tt1

        Tt11 = repmat(t2', size(t1, 1), 1);
        Tt22 = repmat(t1, 1, size(t2, 1));
        Z1 = (Tt11-Tt22)/sigma - sigma*gamma/2;
        z2 = t1/sigma + sigma*gamma/2;
        Z2 = repmat(z2,  1, size(t2, 1));
        [preFactorZ1, gradPreFactorZ1] = quadrantWofzGrad(Z1, sigma,1);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma,1);
        Tt1minusTt2 = Tt11 - Tt22;
        preFactor1 = exp(sigma2*(gamma^2)/4 - gamma*(Tt1minusTt2)...
            + log(sigma2*gamma - 2*(Tt1minusTt2)) + log(gradThetaGamma));
        preFactor2 = exp(-(Tt1minusTt2).*(Tt1minusTt2)/sigma2 ...
            + log(gradThetaGamma) + log(gradPreFactorZ1));
        preFactor3 = exp( repmat(-t1.*t1/sigma2,1, size(t2, 1))...
            - gamma*Tt11 + log(gradThetaGamma) ...
            + log(Tt11.*repmat(preFactorz2,1, size(t2, 1)) + repmat(gradPreFactorz2,1, size(t2, 1))));


    case {3,4}

        % In this mode Tt1 is actually Tt1 and Tt2 is zero.
        % Z1 is a vector and Z2 is a scalar
        z1 = t1/sigma - sigma*gamma/2;
        z2 = sigma*gamma/2;
        Z1 = repmat(z1, 1, size(t2,1));
        Z2 = repmat(z2, size(t1, 1), size(t2,1));
        [preFactorz1, gradPreFactorz1] = quadrantWofzGrad(z1, sigma,1);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma,1);
        preFactor1 = repmat(exp(sigma2*(gamma^2)/4 - gamma*(t1) ...
            + log(sigma2*gamma - 2*t1) + log(gradThetaGamma)), 1, size(t2, 1));
        preFactor2 = repmat(exp(-t1.*t1/sigma2 + log(gradThetaGamma)...
            + log(gradPreFactorz1)), 1, size(t2, 1));
        preFactor3 = repmat(exp(- gamma*t1 + log(gradThetaGamma) ...
            + log(t1*preFactorz2 + gradPreFactorz2)),1, size(t2, 1));

end

g = zeros(size(preFactor1));
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)>=0
ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = preFactor1(ind) - preFactor2(ind) + preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)>=0
ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = preFactor2(ind) + preFactor3(ind);
end
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)<0
ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - preFactor2(ind) - preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)<0
ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = -preFactor1(ind) + preFactor2(ind) - preFactor3(ind);
end;

if mode==4
    g = g.';
end


