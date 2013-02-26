function g = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, pwz1, pwz2, n, m, covDiag)

% SHEATKERNDIAGGRADIENT Gradient of the parameters of diagonal of a SHEAT kernel.
%
%	Description:
%
%	G = SHEATKERNDIAGGRADIENT(SIGMAX, LENGTHX, S, W, GAMMA, PWZ1, PWZ2,
%	N, M, COVDIAG) computes the gradient of the diagonal of a SHEAT
%	kernel parameter.
%	 Returns:
%	  G - gradient of the parameters.
%	 Arguments:
%	  SIGMAX - length-scale of the spatial gp prior.
%	  LENGTHX - length of the spatial domain
%	  S - inputs for which kernel is to be computed.
%	  W - precomputed constant.
%	  GAMMA - precomputed constant.
%	  PWZ1 - precomputed constant.
%	  PWZ2 - precomputed constant.
%	  N - integer indicating first series term
%	  M - integer indicating second series term
%	  COVDIAG - partial derivatives
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, SHEATKERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez


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
