function k = lfmvKernDiagCompute(lfmKern, t)

% LFMVKERNDIAGCOMPUTE Compute diagonal of LFMV kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for velocity- velocity in
% the switching dynamical latent force kernel.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Get length scale out.
sigma2 = 2/lfmKern.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha(1) = lfmKern.damper./(2*lfmKern.mass);
omega(1) = sqrt(lfmKern.spring./lfmKern.mass - alpha(1)*alpha(1));

% Precomputations to increase the speed
preExp1 = zeros(length(t),2);
gamma1_p = alpha(1) + j*omega(1);
gamma1_m = alpha(1) - j*omega(1);
preGamma(1) = gamma1_p + gamma1_p;
preGamma(2) = gamma1_p + gamma1_m;
preGamma(3) = gamma1_m + gamma1_p;
preGamma(4) = gamma1_m + gamma1_m;
preConst = 1./preGamma;
preFactors(1) = preConst(2) - preConst(1);
preFactors(2) = preConst(3) - preConst(4);
preFactors(3) = preConst(3) - preConst(1);
preFactors(4) = preConst(2) - preConst(4);
preExp1(:,1) = gamma1_p.*exp(-gamma1_p*t);
preExp1(:,2) = gamma1_m.*exp(-gamma1_m*t);
% Actual computation of the kernel
sk =lfmDiagComputeH3VV(gamma1_p, gamma1_m, sigma2, t, preFactors([1 2]), 1) + ...
    lfmDiagComputeH3VV(gamma1_p, gamma1_m, sigma2, t, preFactors([3 4]), 0) + ...
    lfmDiagComputeH4VV(gamma1_p, gamma1_m, sigma2, t, preGamma([1 2 4 3]), preExp1 ) + ...
    lfmDiagComputeH4VV(gamma1_p, gamma1_m, sigma2, t, preGamma([1 3 4 2]), preExp1 );

if lfmKern.isNormalised
    k0 = lfmKern.sensitivity^2/(8*sqrt(2)*lfmKern.mass^2*omega^2);
else
    k0 = sqrt(pi)*sigma*lfmKern.sensitivity^2/(8*lfmKern.mass^2*omega^2);
end
k = k0*sk;


function h = lfmDiagComputeH3VV(gamma1_p, gamma1_m, sigma2, t, preFactor, mode)

h = preFactor(1)*lfmvvComputeUpsilonDiagVector(gamma1_p, sigma2, t, mode) ...
    + preFactor(2)*lfmvvComputeUpsilonDiagVector(gamma1_m,sigma2, t, mode);

function h = lfmDiagComputeH4VV(gamma1_p, gamma1_m, sigma2, t, ...
    preFactor, preExp)

h =   lfmvpComputeUpsilonVector(gamma1_p,sigma2, t, 0).*( preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)) ...
    + lfmvpComputeUpsilonVector(gamma1_m,sigma2, t, 0).*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3));




