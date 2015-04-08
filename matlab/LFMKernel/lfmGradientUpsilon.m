function    g = lfmGradientUpsilon(gamma,sigma2,gradThetaGamma,Tt1,Tt2);

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
% RETURN g : Gradient of the kernel with respect to the desired
% parameter.
%
% COPYRIGHT : David Luengo, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientH

% LFM


%Parameters of the function

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

g = zeros(size(Tt1));

Z1 = (Tt1-Tt2)/sigma - sigma*gamma/2;
Z2 = Tt2/sigma + sigma*gamma/2;

%%% Gradient Evaluation %%%

% Evaluation of the gradient when real(Z1)>=0 and real(Z2)>=0

ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = exp(sigma2*(gamma^2)/4 - gamma*(Tt1(ind)-Tt2(ind)) ...
            + log(sigma2*gamma - 2*(Tt1(ind)-Tt2(ind))) + log(gradThetaGamma)) ...
        - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(sigma*gradThetaGamma) + log(1/sqrt(pi) - Z1(ind).*W(j*Z1(ind)))) ...
        + exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(gradThetaGamma) ...
            + log(Tt1(ind).*W(j*Z2(ind)) + sigma*(1/sqrt(pi) - Z2(ind).*W(j*Z2(ind)))));
end

% Evaluation of the gradient when real(Z1)<0 and real(Z2)>=0

ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(sigma*gradThetaGamma) + log(1/sqrt(pi) + Z1(ind).*W(-j*Z1(ind)))) ...
        + exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(gradThetaGamma) ...
            + log(Tt1(ind).*W(j*Z2(ind)) + sigma*(1/sqrt(pi) - Z2(ind).*W(j*Z2(ind)))));
end

% Evaluation of the gradient when real(Z1)>=0 and real(Z2)<0

ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(sigma*gradThetaGamma) + log(1/sqrt(pi) - Z1(ind).*W(j*Z1(ind)))) ...
        - exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(gradThetaGamma) ...
            + log(Tt1(ind).*W(-j*Z2(ind)) - sigma*(1/sqrt(pi) + Z2(ind).*W(-j*Z2(ind)))));
end

% Evaluation of the gradient when real(Z1)<0 and real(Z2)<0

ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - exp(sigma2*(gamma^2)/4 - gamma*(Tt1(ind)-Tt2(ind)) ...
            + log(sigma2*gamma - 2*(Tt1(ind)-Tt2(ind))) + log(gradThetaGamma)) ...
        - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(sigma*gradThetaGamma) + log(1/sqrt(pi) + Z1(ind).*W(-j*Z1(ind)))) ...
        + exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(gradThetaGamma) ...
            + log(Tt1(ind).*W(j*Z2(ind)) + sigma*(1/sqrt(pi) - Z2(ind).*W(j*Z2(ind)))));
end
