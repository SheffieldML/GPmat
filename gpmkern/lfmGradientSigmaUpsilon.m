function g = lfmGradientSigmaUpsilon(gamma, sigma2, t1, t2, mode)

% LFMGRADIENTSIGMAUPSILON Gradient of the function \upsilon(z) with respect
% to \sigma.
% FORMAT
% DESC Computes the gradient of the function \upsilon(z) with respect to
% the length-scale of the input "force", \sigma.
% ARG gamma : Gamma value of the system.
% ARG sigma2 : length scale of latent process.
% ARG Tt1 : first time input (number of time points 1 x number of time points 2).
% ARG Tt2 : second time input (number of time points 1 x number of time points 2).
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN g : Gradient of the kernel with respect to the desired
% parameter.
%
% COPYRIGHT : David Luengo, 2008
%
% MODIFICATIONS: Mauricio Alvarez, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaH

% KERN


%Parameters of the function

sigma = sqrt(sigma2);

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
        [preFactorZ1, gradPreFactorZ1] = quadrantWofzGrad(Z1, sigma, 0);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma, 0);
        Tt1minusTt2 = Tt11 - Tt22;
        preFactor1 = sigma*(gamma^2)*exp(sigma2*(gamma^2)/4 - gamma*(Tt1minusTt2));
        preFactor2 = 2*exp(-(Tt1minusTt2).*(Tt1minusTt2)/sigma2 + log((Tt1minusTt2).*(Tt1minusTt2).*preFactorZ1/(sigma^3) ...
            + ((Tt1minusTt2)/sigma2+gamma/2).* (gradPreFactorZ1)));
        preFactor3 = 2*exp(repmat((-t2.*t2/sigma2 + log(t2.*t2.*preFactorz2/(sigma^3) ...
            + (t2/sigma2-gamma/2).* gradPreFactorz2)).',size(t1, 1), 1)-gamma*Tt11);

    case 2
        % In this mode the Tt1 is actually the Tt2 and Tt2 is actually Tt1
        % Z1 is a matrix with the order given by Tt2 and z2 is a vector
        % with the order given by Tt1
        Tt11 = repmat(t2', size(t1, 1), 1);
        Tt22 = repmat(t1, 1, size(t2, 1));
        Z1 = (Tt11-Tt22)/sigma - sigma*gamma/2;
        z2 = t1/sigma + sigma*gamma/2;
        Z2 = repmat(z2,  1, size(t2, 1));
        [preFactorZ1, gradPreFactorZ1] = quadrantWofzGrad(Z1, sigma, 0);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma, 0);
        Tt1minusTt2 = Tt11 - Tt22;
        preFactor1 = sigma*(gamma^2)*exp(sigma2*(gamma^2)/4 - gamma*(Tt1minusTt2));
        preFactor2 = 2*exp(-(Tt1minusTt2).*(Tt1minusTt2)/sigma2 + log((Tt1minusTt2).*(Tt1minusTt2).*preFactorZ1/(sigma^3) ...
            + ((Tt1minusTt2)/sigma2+gamma/2).* (gradPreFactorZ1)));
        preFactor3 = 2*exp(repmat((-t1.*t1/sigma2 + log(t1.*t1.*preFactorz2/(sigma^3) ...
            + (t1/sigma2-gamma/2).* gradPreFactorz2)),1, size(t2, 1))-gamma*Tt11);
    case {3,4}
        % In this mode Tt1 is actually Tt1 and Tt2 is zero.
        % Z1 is a vector and Z2 is a scalar
        z1 = t1/sigma - sigma*gamma/2;
        z2 = sigma*gamma/2;
        Z1 = repmat(z1, 1, size(t2,1));
        Z2 = repmat(z2, size(t1, 1), size(t2,1));
        [preFactorz1, gradPreFactorz1] = quadrantWofzGrad(z1, sigma,0);
        [preFactorz2, gradPreFactorz2] = quadrantWofzGrad(z2, sigma,0);
        preFactor1 = repmat(sigma*(gamma^2)*exp(sigma2*(gamma^2)/4 - gamma*t1), 1, size(t2, 1));
        preFactor2 = repmat(2*exp(-t1.*t1/sigma2 + log(t1.*t1.*preFactorz1/(sigma^3) ...
            + ((t1)/sigma2+gamma/2).*(gradPreFactorz1))),1, size(t2, 1));
        preFactor3 = 2*repmat( exp( log(-gamma/2.* gradPreFactorz2)-gamma*t1) ,1,size(t2, 1));
end

g = zeros(size(preFactor1));
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)>=0
ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = preFactor1(ind) - preFactor2(ind) - preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)>=0
ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = preFactor2(ind) - preFactor3(ind);
end
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)<0
ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - preFactor2(ind) + preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)<0
ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = -preFactor1(ind) + preFactor2(ind) + preFactor3(ind);
end;

if mode==4
    g = g.';
end
