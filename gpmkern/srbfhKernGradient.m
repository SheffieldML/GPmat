function g = srbfhKernGradient(sigmax, lengthX, s1, s2, w, gamma, n, covGrad, wz1, wz2)

% SRBFHKERNGRADIENT Gradient of the parameters of a SRBFH kernel.
% FORMAT
% DESC computes the gradient of a SRBFH kernel parameter.
% ARG sigmax : length-scale of the spatial gp prior.
% ARG lengthX : length of the spatial domain
% ARG s1 : row inputs for which kernel is to be computed. 
% ARG s2 : column inputs for which kernel is to be computed. 
% ARG w  : precomputed constant.
% ARG gamma : precomputed constant.
% ARG n : integer indicating first series term
% ARG covGrad : partial derivatives
% ARG wz1 : precomputed factors
% ARG wz2 : precomputed factors
% RETURN g : gradient of the parameters.
%
% SEEALSO : multiKernParamInit, multiKernCompute, srbfhKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sinS1 = sin(w(n)*s1);
bTerm = sigmax*gamma(n)/2;
argz2 = (s2-lengthX)/sigmax;
argz1 =  s2/sigmax;

if nargin < 9
    z1 = argz1 + bTerm;
    z2 = argz2 + bTerm;
    wz1 = wofzPoppe(sqrt(-1)*z1);
    wz2 = wofzPoppe(sqrt(-1)*z2);
end

dWdsx = (2*((s2 - lengthX).^2)/(sigmax^3) + (2/sigmax)*(bTerm^2  - argz2.^2))...
    .*exp(-argz2.^2 + gamma(n)*lengthX).*wz2 - (2/(sigmax*sqrt(pi)))...
    *(bTerm-argz2).*exp(-argz2.^2 + gamma(n)*lengthX) ...
    - (2*(s2.^2)/(sigmax^3) + (2/sigmax)*(bTerm^2  - argz1.^2)).*exp(-argz1.^2).*wz1 ...
    + (2/(sigmax*sqrt(pi)))*(bTerm-argz1).*exp(-argz1.^2);

W = exp(-argz2.^2 + gamma(n)*lengthX).*wz2 -  exp(-argz1.^2).*wz1;

dvec = (sqrt(pi)/2)*imag(W) + (sigmax*sqrt(pi)/2)*imag(dWdsx);


g = sum(sum((sinS1*dvec').*covGrad));


