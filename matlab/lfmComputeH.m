function  h = lfmComputeH(gamma1, gamma2, sigma2, t1, t2, mode)

% LFMCOMPUTEH Helper function for computing part of the LFM kernel.
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
% COPYRIGHT : David Luengo, 2007, Mauricio Alvarez, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007
%
%
%
% SEEALSO : lfmKernParamInit, lfmXlfmKernCompute

% KERN

% Evaluation of h


if mode==1
    % t1 is really t1 and t2 is really t2
%    Tt1 = repmat(t2', size(t1, 1), 1);
    Tt1 = repmat(t1, 1, size(t2, 1));
    h = (lfmComputeUpsilon(gamma1,sigma2,t2,t1,1) - exp(-gamma2*Tt1) ...
        .* lfmComputeUpsilon(gamma1,sigma2,t2, zeros(size(t1)),3))/(gamma1+gamma2);
else
    % t1 is really t2 and t2 is really t1
    Tt1 = repmat(t1', size(t2, 1), 1);
    %Tt1 = repmat(t2', size(t1, 1), 1);
%    Tt1 = repmat(t1, 1, size(t2, 1))';
    h = (lfmComputeUpsilon(gamma1,sigma2,t2,t1,2) - exp(-gamma2*Tt1) ...
        .* lfmComputeUpsilon(gamma1,sigma2,t2, zeros(size(t1)),4))/(gamma1+gamma2);
end

