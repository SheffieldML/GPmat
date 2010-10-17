function g = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, pwz1, pwz2, n, m, covGrad)

% SHEATKERNGRADIENT Gradient of the parameters of a SHEAT kernel.
% FORMAT
% DESC computes the gradient of a SHEAT kernel parameter.
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
% ARG covGrad : partial derivatives
% RETURN g : gradient of the parameters.
%
% SEEALSO : multiKernParamInit, multiKernCompute, sheatKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sinS1 = sin(w(n)*s1);
sinS2 = sin(w(m)*s2);
if n == m
    wn = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
    dwn = ((sigmax*gamma(n)^2)/2)*pwz1(n) - gamma(n)/(sqrt(pi)) ...
        - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n)/2)^2 - (lengthX/sigmax)^2))...
        *exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n) + ...
        (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
        *exp(-gamma(n)*lengthX);
    dC = (lengthX*sqrt(pi)/2)*real(wn) + (sigmax*lengthX*sqrt(pi)/2)*real(dwn);
    if mod(n,2) == 0
        dC = dC + sigmax*(exp(-(lengthX/sigmax)^2) - 1) + ((lengthX^2)/sigmax)*exp(-(lengthX/sigmax)^2);
    else
        dC = dC - sigmax*(exp(-(lengthX/sigmax)^2) + 1) - ((lengthX^2)/sigmax)*exp(-(lengthX/sigmax)^2);
    end
else
    if mod(n+m,2)==1
        wm = pwz1(m) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m);
        wn = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
        dwn = ((sigmax*gamma(n)^2)/2)*pwz1(n) - gamma(n)/(sqrt(pi)) ...
            - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n)/2)^2 - (lengthX/sigmax)^2))...
            *exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n) + ...
            (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            *exp(-gamma(n)*lengthX);
        dwm = ((sigmax*gamma(m)^2)/2)*pwz1(m) - gamma(m)/(sqrt(pi)) ...
            - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(m)/2)^2 - (lengthX/sigmax)^2))...
            *exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m) + ...
            (2/(sigmax*sqrt(pi)))*(sigmax*gamma(m)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            *exp(-gamma(m)*lengthX);
        A = (1/(n+m)) + (1/(n-m));
        B = (1/(n+m)) - (1/(n-m));
        dC = (lengthX/(2*sqrt(pi)))*(A*imag(wm) + B*imag(wn)) ...
            + (sigmax*lengthX/(2*sqrt(pi)))*(A*imag(dwm) + B*imag(dwn));
     else
        dC = 0;
    end
end

g = dC*sum(sum((sinS1*sinS2').*covGrad));

