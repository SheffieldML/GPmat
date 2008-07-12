function upsilon = lfmComputeUpsilon(gamma,sigma2,Tt1,Tt2);

% LFMCOMPUTEUPSILON Helper function for comptuing part of the LFM kernel.
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma : Gamma value for the system.
% ARG sigma2 : length scale of latent process.
% ARG Tt1 : first time input (number of time points of x1 x number of time
% points of x2).
% ARG Tt2 : second time input (number of time points x number of time
% points of x2).
% RETURN upsilon : result of this subcomponent of the kernel for the given
% values.
%
%
% COPYRIGHT : David Luengo, 2008
%
% MODIFICATIONS  : Mauricio Alvarez, 2008
%
% SEEALSO : lfmKernParamInit, lfmXlfmKernCompute, lfmComputeH, W

% KERN


dev = sqrt(sigma2);

% Initialization of vectors and matrices

upsilon = zeros(size(Tt1));

Z1 = (Tt1-Tt2)/dev - dev*gamma/2;
Z2 = Tt2/dev + dev*gamma/2;


%%% Evaluation of Upsilon %%%

% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)>=0

ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    upsilon(ind) = 2*exp(sigma2*(gamma^2)/4 - gamma*(Tt1(ind)-Tt2(ind))) ...
        - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 + log(wofz(j*Z1(ind)))) ...
        - exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(wofz(j*Z2(ind))));
end;

% Evaluation of Upsilon when real(Z1)<0 and real(Z2)>=0

ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    upsilon(ind) = exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(wofz(-j*Z1(ind)))) ...
        - exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) + log(wofz(j*Z2(ind))));
end;

% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)<0

ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    upsilon(ind) = - exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(wofz(j*Z1(ind)))) ...
        + exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) ...
            + log(wofz(-j*Z2(ind))));
end;

% Evaluation of Upsilon when real(Z1)<0 and real(Z2)<0

ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    upsilon(ind) = exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log(wofz(-j*Z1(ind)))) ...
        - 2*exp(sigma2*(gamma^2)/4-gamma*(Tt1(ind)-Tt2(ind))) ...
        + exp(-Tt2(ind).*Tt2(ind)/sigma2 - gamma*Tt1(ind) ...
            + log(wofz(-j*Z2(ind))));
end;

return;
