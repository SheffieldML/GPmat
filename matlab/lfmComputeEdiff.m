function x = lfmComputeEdiff(gamma,sigma2,Tt,Tt2);
%
%
% Function that evaluates the function E_i, used in Kxf and Kxx. Requires
% the function erfz, which provides the erf(z) for complex values.
%
%   Hi = GPLfmEvalEi(gamma,sigma2,Tt,Tt2);
%
%   Ei     - Matrix TxT2 with the desired values of E_i.
%   gamma  - Propagation constant.
%   sigma2 - Length scale of the GP input.
%   Tt     - Matrix TxT2. Time instants for x_k.
%   Tt2    - Matrix TxT2. Time instants for x_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 29 October 2007
% Last Modification : 29 October 2007


% Parameters of the kernel

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

x = zeros(size(Tt));

% Evaluation of Ei

x = erfz((Tt-Tt2)/sigma-sigma*gamma/2) + erfz(Tt2/sigma+sigma*gamma/2);
