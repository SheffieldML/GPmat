function  [h, compUpV, compUp] =  lfmComputeH4VP(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp, mode)

% LFMCOMPUTEH4VP Helper function for computing part of the LFMVXLFM kernel.
% FORMAT
% DESC computes a portion of the LFMVXLFM kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO : lfmComputeH4Hat, lfmXlfmKernCompute

% KERN

% This could also be used with str2func changing between 'lfm' and 'lfmvp'


if mode==0
    if nargout > 1
        [compUpV{1}, compUp{1}] = lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1, 0);
        [compUpV{2}, compUp{2}] = lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1, 0);
        h =  compUpV{1}*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + compUpV{2}*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    else
        h =  lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1, 0)*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1, 0)*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    end
else
    if nargout > 1
        compUp{1} = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUp{2} = lfmComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUp{1}*(  preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + compUp{2}*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
        compUpV = compUp;
    else
        h =  lfmComputeUpsilonVector(gamma1_p,sigma2, t1)*(  preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + lfmComputeUpsilonVector(gamma1_m,sigma2, t1)*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    end
end
