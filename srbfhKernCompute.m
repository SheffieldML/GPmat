function [K, wz1, wz2] = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, n)

% SRBFHKERNCOMPUTE Compute an SRBFH kernel.
% FORMAT

%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sinS1 = sin(w(n)*s1);
bTerm = sigmax*gamma(n)/2;
argz2 = (s2-lengthX)/sigmax;
z2 = argz2 + bTerm;
argz1 =  s2/sigmax;
z1 = argz1 + bTerm;
wz2 = wofzPoppe(sqrt(-1)*z2);
wz1 = wofzPoppe(sqrt(-1)*z1);
vecXp = (sigmax*sqrt(pi)/2)*imag(exp(-argz2.^2 + gamma(n)*lengthX).*wz2 -  exp(-argz1.^2).*wz1);
K = sinS1*vecXp';

