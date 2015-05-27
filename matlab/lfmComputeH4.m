function  [h, compUp] =  lfmComputeH4(gamma1_p, gamma1_m, sigma2, t1, preFactor, preExp,...
    mode, term )
% LFMCOMPUTEH4 Helper function for computing part of the LFM kernel.
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : David Luengo, 2007
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007
%
%
%
% SEEALSO : lfmKernParamInit, lfmXlfmKernCompute

% KERN

% Evaluation of h

if nargin<8
    term =[];
end

if ~mode
    if ~term
        if nargout>1
            compUp = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
            h = compUp*( preExp/preFactor(1) - conj(preExp)/preFactor(2)).';
        else
            h = lfmComputeUpsilonVector(gamma1_p,sigma2, t1)*( preExp/preFactor(1) - conj(preExp)/preFactor(2)).';
        end
    else
        if nargout>1
            compUp = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
            h = compUp*(preExp/preFactor(1)).' - conj(compUp)*(preExp/preFactor(2)).';
        else
            upsilon = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
            h = upsilon*(preExp/preFactor(1)).' - conj(upsilon)*(preExp/preFactor(2)).';
        end
    end
else
    if nargout > 1
        compUp{1} = lfmComputeUpsilonVector(gamma1_p,sigma2, t1);
        compUp{2} = lfmComputeUpsilonVector(gamma1_m,sigma2, t1);
        h =  compUp{1}*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + compUp{2}*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    else
        h =  lfmComputeUpsilonVector(gamma1_p,sigma2, t1)*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
            + lfmComputeUpsilonVector(gamma1_m,sigma2, t1)*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
    end
end
