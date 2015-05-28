function g = srbfhKernGradient(sigmax, lengthX, s1, s2, w, gamma, n, covGrad, wz1, wz2)

% SRBFHKERNGRADIENT Gradient of the parameters of a SRBFH kernel.
%
%	Description:
%
%	G = SRBFHKERNGRADIENT(SIGMAX, LENGTHX, S1, S2, W, GAMMA, N, COVGRAD,
%	WZ1, WZ2) computes the gradient of a SRBFH kernel parameter.
%	 Returns:
%	  G - gradient of the parameters.
%	 Arguments:
%	  SIGMAX - length-scale of the spatial gp prior.
%	  LENGTHX - length of the spatial domain
%	  S1 - row inputs for which kernel is to be computed.
%	  S2 - column inputs for which kernel is to be computed.
%	  W - precomputed constant.
%	  GAMMA - precomputed constant.
%	  N - integer indicating first series term
%	  COVGRAD - partial derivatives
%	  WZ1 - precomputed factors
%	  WZ2 - precomputed factors
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, SRBFHKERNCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez


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


