function  [h, compUp] =  lfmComputeH4JA(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp, mode)

% LFMCOMPUTEH4JA Helper function for computing part of the LFMJA kernel.
% FORMAT
% DESC computes a portion of the LFMJA kernel.
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
% SEEALSO : lfmComputeH4.m, lfmComputeH4AA.m

% KERN

if mode==0
    if nargout > 1
        compUp{1} = lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUp{2} = lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUp{1}*(preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + compUp{2}*(preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    else
        h =  lfmjpComputeUpsilonVector(gamma1_p,sigma2, t1)*(preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + lfmjpComputeUpsilonVector(gamma1_m,sigma2, t1)*(preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    end
else
    if nargout > 1
        compUp{1} = lfmapComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUp{2} = lfmapComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUp{1}*(preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + compUp{2}*(preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    else
        h =  lfmapComputeUpsilonVector(gamma1_p,sigma2, t1)*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
            + lfmapComputeUpsilonVector(gamma1_m,sigma2, t1)*(preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
    end
end
