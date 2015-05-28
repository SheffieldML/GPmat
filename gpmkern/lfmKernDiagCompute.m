function k = lfmKernDiagCompute(kern, t)

% LFMKERNDIAGCOMPUTE Compute diagonal of LFM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the single input
% motif kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : lfmKernParamInit, kernDiagCompute, kernCreate, lfmKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% COPYRIGHT : Mauricio A. Alvarez, 2008, 2010

% KERN

if size(t, 2) > 1 
  error('Input can only have one column');
end

% Get length scale out.
sigma2 = 2/kern.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha = kern.damper./(2*kern.mass);
omega = sqrt(kern.spring./kern.mass - alpha*alpha);

if isreal(omega)
    % Precomputations to increase speed
    gamma  = alpha + j*omega;
    preExp =  exp(-gamma*t);
    % Actual computation of the kernel
    sk = real(lfmDiagComputeH3(-gamma, sigma2, t, [2*alpha gamma], preExp, 0) + ...
        lfmDiagComputeH4( gamma, sigma2, t, [gamma 2*alpha], preExp ,0));
    if  kern.isNormalised
        k0 = kern.sensitivity^2/(4*sqrt(2)*kern.mass^2*omega^2);
    else
        k0 = sqrt(pi)*sigma*kern.sensitivity^2/(4*kern.mass^2*omega^2);
    end
    k = k0*sk;  
else
    % Precomputations to increase the speed
    preExp = zeros(length(t),2);    
    gamma_p = alpha + j*omega;
    gamma_m = alpha - j*omega;
    preFactors(1) = 2/(gamma_p + gamma_m) - 1/gamma_m;
    preFactors(2) = 2/(gamma_p + gamma_m) - 1/gamma_p;
    preExp(:,1) = exp(-gamma_p*t);
    preExp(:,2) = exp(-gamma_m*t);
    % Actual computation of the kernel    
    sk = lfmDiagComputeH3(-gamma_m, sigma2, t, preFactors(1), preExp(:,2), 1) + ...
         lfmDiagComputeH3(-gamma_p, sigma2, t, preFactors(2), preExp(:,1), 1) + ...
         lfmDiagComputeH4( gamma_m, sigma2, t, [gamma_m  (gamma_p + gamma_m)], [preExp(:,2) preExp(:,1)] , 1) + ...
         lfmDiagComputeH4( gamma_p, sigma2, t, [gamma_p  (gamma_p + gamma_m)], preExp , 1);

    if kern.isNormalised
        k0 = kern.sensitivity^2/(8*sqrt(2)*kern.mass^2*omega^2);
    else
        k0 = sqrt(pi)*sigma*kern.sensitivity^2/(8*kern.mass^2*omega^2);
    end
    k = k0*sk; 
    
end
end

function vec = lfmDiagComputeH3(gamma, sigma2, t, factor, preExp, mode)

if mode
    vec = preExp.*lfmComputeUpsilonVector(gamma ,sigma2, t)*factor;
else
    temp = preExp.*lfmComputeUpsilonVector(gamma ,sigma2, t);    
    vec = 2*real(temp/factor(1)) - temp/factor(2);
end
end

function vec = lfmDiagComputeH4(gamma, sigma2, t, factor, preExp, mode)

if mode
    vec = (preExp(:,1)/factor(1) -  2*preExp(:,2)/factor(2))...
        .*lfmComputeUpsilonVector(gamma ,sigma2, t);
else
    temp = lfmComputeUpsilonVector(gamma ,sigma2, t);
    temp2 = temp.*conj(preExp)/factor(2);
    vec = temp.*preExp/factor(1) - 2*real(temp2);
end

end




