function g = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, pwz1, pwz2, n, m, covDiag)

% SHEATKERNDIAGGRADIENT Gradient of the parameters of diagonal of a SHEAT kernel.
% FORMAT
% DESC computes the gradient of the diagonal of a SHEAT kernel parameter.
% ARG sigmax : length-scale of the spatial gp prior.
% ARG lengthX : length of the spatial domain
% ARG s : inputs for which kernel is to be computed.
% ARG w : precomputed constant.
% ARG gamma : precomputed constant.
% ARG pwz1 : precomputed constant.
% ARG pwz2 : precomputed constant.
% ARG n : integer indicating first series term
% ARG m : integer indicating second series term
% ARG covDiag : partial derivatives
% RETURN g : gradient of the parameters.
%
% SEEALSO : multiKernParamInit, multiKernCompute, sheatKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sinS1 = sin(w(n)*s);
sinS2 = sin(w(m)*s);

if n == m
    Wox = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
    
    dWox = ((sigmax*gamma(n)^2)/2)*pwz1(n) - gamma(n)/(sqrt(pi)) ...
        - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n)/2)^2 - (lengthX/sigmax)^2))...
        *exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n) + ...
        (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
        *exp(-gamma(n)*lengthX);
    
    dC =  (sqrt(pi)*lengthX/2)*(real(Wox) ...
        - imag(Wox)*((sigmax^2*n*pi)/(2*lengthX^2) + (1/(n*pi)))) ...
        + (sigmax*sqrt(pi)*lengthX/2)*(real(dWox) ...
        - imag(dWox)*((sigmax^2*n*pi)/(2*lengthX^2) + (1/(n*pi))) ...   
        - imag(Wox)*((sigmax*n*pi)/(lengthX^2))) ...
        + sigmax*(exp(-(lengthX/sigmax)^2)*cos(n*pi) - 1)  ...
        + ((lengthX^2)/sigmax)*exp(-(lengthX/sigmax)^2)*cos(n*pi);
else
    if mod(n+m,2)==1
        dC = 0;
    else
        Woxm = pwz1(m) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m);
        Woxn = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
        
        dWoxm = ((sigmax*gamma(m)^2)/2)*pwz1(m) - gamma(m)/(sqrt(pi)) ...
            - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(m)/2)^2 - (lengthX/sigmax)^2))...
            *exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m) + ...
            (2/(sigmax*sqrt(pi)))*(sigmax*gamma(m)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            *exp(-gamma(m)*lengthX);
        dWoxn = ((sigmax*gamma(n)^2)/2)*pwz1(n) - gamma(n)/(sqrt(pi)) ...
            - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n)/2)^2 - (lengthX/sigmax)^2))...
            *exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n) + ...
            (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            *exp(-gamma(n)*lengthX);
                
        dC = (lengthX/(sqrt(pi)*(m^2-n^2)))*(n*imag(Woxm) - m*imag(Woxn)) ...
            + ((sigmax*lengthX)/(sqrt(pi)*(m^2-n^2)))*(n*imag(dWoxm) - m*imag(dWoxn));      
    end
end

g = dC*sum((sinS1.*sinS2).*covDiag);
