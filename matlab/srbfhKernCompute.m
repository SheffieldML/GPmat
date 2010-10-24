function [K, wz1, wz2] = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, n)

% SHEATKERNCOMPUTE Compute a cross kernel between two SHEAT kernels.
% FORMAT
% DESC computes cross kernel terms between two SHEAT kernels. It gives a 
% partial result needed for the computation of the HEAT kernel. 
% ARG sigmax : length-scale of the spatial gp prior.
% ARG lengthX : length of the spatial domain
% ARG s1 : row inputs for which kernel is to be computed. 
% ARG s2 : column inputs for which kernel is to be computed. 
% ARG w  : precomputed constant.
% ARG gamma : precomputed constant.
% ARG n : integer indicating first series term
% RETURN K : block of values from kernel matrix.
% RETURN wz1 : precomputed values of the complex error function.
% RETURN wz2 : precomputed values of the complex error function.
%
% SEEALSO : multiKernParamInit, multiKernCompute, heatKernParamInit
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

