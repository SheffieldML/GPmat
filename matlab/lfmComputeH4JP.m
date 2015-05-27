function  [h, compUpJP] =  lfmComputeH4JP(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp, mode)

% LFMCOMPUTEH4JP Helper function for computing part of the LFMJP kernel.
% FORMAT
% DESC computes a portion of the LFMAP kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmComputeH4.m, lfmComputeH4AP.m

% KERN

if mode==0
    if nargout > 1
        compUpJP{1} = lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUpJP{2} = lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUpJP{1}*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + compUpJP{2}*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    else
        h =  lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1)*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1)*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    end
else
    if nargout > 1
        compUp{1} = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUp{2} = lfmComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUp{1}*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + compUp{2}*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
        compUpJP = compUp;
    else
        h =  lfmComputeUpsilonVector(gamma1_p,sigma2, t1)*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1) ).' ...
            + lfmComputeUpsilonVector(gamma1_m,sigma2, t1)*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    end
end
