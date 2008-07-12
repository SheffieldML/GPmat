function  h = lfmComputeH(gamma1, gamma2, sigma2, Tt1, Tt2)

% LFMCOMPUTEH Helper function for computing part of the LFM kernel.
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : David Luengo, 2007, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007
%
% MODIFICATIONS : Mauricio Alvarez, 2008
%
%
% SEEALSO : lfmKernParamInit, lfmXlfmKernCompute

% KERN

% Evaluation of h

h = (lfmComputeUpsilon(gamma1,sigma2,Tt2,Tt1) - exp(-gamma2*Tt1) ...
    .* lfmComputeUpsilon(gamma1,sigma2,Tt2,zeros(size(Tt1))))/(gamma1+gamma2);

