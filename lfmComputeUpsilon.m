function upsilon = lfmComputeUpsilon(gamma,sigma2,t1,t2, mode)

% LFMCOMPUTEUPSILON Helper function for comptuing part of the LFM kernel.
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma : Gamma value for the system.
% ARG sigma2 : length scale of latent process.
% ARG Tt1 : first time input (number of time points of x1 x number of time
% points of x2).
% ARG Tt2 : second time input (number of time points x number of time
% points of x2).
% ARG mode : indicates in which way the variables t1 and t2 must be
% transposed
% RETURN upsilon : result of this subcomponent of the kernel for the given
% values.
%
%
% COPYRIGHT : David Luengo, Mauricio Alvarez, 2008
%
%
% SEEALSO : lfmKernParamInit, lfmXlfmKernCompute, lfmComputeH, W

% KERN

dev = sqrt(sigma2);

switch mode
    case {1, 2}
        % In this mode the Tt1 is actually the Tt1 and Tt2 is actually Tt2
        % Z1 is a matrix with order given by Tt1 and z2 is a vector with
        % order given by Tt2
        Tt11 = repmat(t1, 1, size(t2, 1));
        Tt22 = repmat(t2', size(t1, 1), 1);
        Z1 = (Tt11-Tt22)/dev - dev*gamma/2;
        z2 = t2/dev + dev*gamma/2;
        Z2 = repmat(z2', size(t1, 1), 1);
%         preFactorZ1 = log(quadrantWofz(Z1));
%         preFactorz2 = log(quadrantWofz(z2));
        preFactorZ1 = logWofzHui(j*Z1);
        preFactorz2 = logWofzHui(j*z2);
        Tt1minusTt2 = Tt11 - Tt22;
        preFactor1 = 2*exp(sigma2*(gamma^2)/4 - gamma*(Tt1minusTt2));
        preFactor2 = exp(-(Tt1minusTt2).*(Tt1minusTt2)/sigma2 + preFactorZ1);
        preFactor3 = exp(repmat((-t2.*t2/sigma2 + preFactorz2).',size(t1, 1), 1) - gamma*Tt11);        
    case {3,4}

        % In this mode Tt1 is actually Tt1 and Tt2 is zero.
        % Z1 is a vector and Z2 is a scalar equal to zero
        z1 = t1/dev - dev*gamma/2;
        z2 = dev*gamma/2;
        Z1 = repmat(z1, 1, size(t2, 1));
        Z2 = repmat(z2, length(t1),length(t2));
%         preFactorz1 = log(quadrantWofz(z1));
%         preFactorz2 = log(quadrantWofz(z2));
        preFactorz1 = logWofzHui(j*z1);
        preFactorz2 = logWofzHui(j*z2);
        preFactor1 = repmat(2*exp(sigma2*(gamma^2)/4 - gamma*(t1)),1, size(t2, 1));
        preFactor2 = repmat(exp(-t1.*t1/sigma2 + preFactorz1),1, size(t2, 1));
        preFactor3 = repmat(exp( - gamma*t1 + preFactorz2),1, size(t2, 1));
end

upsilon = zeros(size(preFactor1));
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)>=0
ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    upsilon(ind) = preFactor1(ind) - preFactor2(ind) - preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)>=0
ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    upsilon(ind) =  preFactor2(ind) - preFactor3(ind);
end
% Evaluation of Upsilon when real(Z1)>=0 and real(Z2)<0
ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    upsilon(ind) = - preFactor2(ind) + preFactor3(ind);
end;
% Evaluation of Upsilon when real(Z1)<0 and real(Z2)<0
ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    upsilon(ind) = -preFactor1(ind) + preFactor2(ind) + preFactor3(ind);
end;

if ((mode == 1) || (mode == 3))
    upsilon = upsilon.';
end



