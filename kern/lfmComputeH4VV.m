function  [h, compUpVP, compUp] =  lfmComputeH4VV(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp)

% LFMCOMPUTEH4VV Helper function for computing part of the LFMVXLFMV kernel.
% FORMAT
% DESC computes a portion of the LFMVXLFMV kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG preFactor : precomputed constants.
% ARG preExp : precomputed exponentials.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO :  lfmvXlfmvKernCompute.m

% KERN


if nargout > 1
    [compUpVP{1}, compUp{1}] = lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1, 0);
    [compUpVP{2}, compUp{2}] = lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1, 0);
    h =  compUpVP{1}*( preExp(:,2)/preFactor(2)- preExp(:,1)/preFactor(1) ).' ...
        + compUpVP{2}*(  preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
else
    h =   lfmvpComputeUpsilonVector(gamma1_p,sigma2, t1, 0)*( preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
        + lfmvpComputeUpsilonVector(gamma1_m,sigma2, t1, 0)*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
end
