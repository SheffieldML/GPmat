function [K, const] = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, pwz1, pwz2, n, m)

% SHEATKERNCOMPUTE Compute a cross kernel between two SHEAT kernels.
% FORMAT
% DESC computes cross kernel terms between two SHEAT kernels. It gives a 
% partial result needed for the computation of the HEAT kernel. 
% ARG sigmax : length-scale of the spatial gp prior.
% ARG lengthX : length of the spatial domain
% ARG s1 : row inputs for which kernel is to be computed. 
% ARG s2 : column inputs for which kernel is to be computed. 
% ARG w : precomputed constant.
% ARG gamma : precomputed constant.
% ARG pwz1 : precomputed constant. 
% ARG pwz2 : precomputed constant. 
% ARG n : integer indicating first series term
% ARG m : integer indicating second series term
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, heatKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sinS1 = sin(w(n)*s1);
sinS2 = sin(w(m)*s2);

if n == m
    const = real(pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n));
    const = (sigmax*lengthX*sqrt(pi)/2)*const; 
    if mod(n,2) == 0
        const = const + (sigmax^2)/2*(exp(-(lengthX/sigmax)^2) - 1);
    else
        const = const - (sigmax^2)/2*(exp(-(lengthX/sigmax)^2) + 1);
    end    
else
    if mod(n+m,2)==1   
        A = (1/(n+m)) + (1/(n-m)); 
        B = (1/(n+m)) - (1/(n-m));
        preA = pwz1(m) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m);
        preB = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
        const = A*imag(preA) + B*imag(preB);
        const = (sigmax*lengthX/(2*sqrt(pi)))*const;
    else
        const = 0;
    end
end

constSinS1 = const*sinS1;
K = constSinS1*sinS2';

