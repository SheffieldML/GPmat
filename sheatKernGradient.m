function g = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, pwz1, pwz2, n, m, covGrad, bctype, isNumeric)

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
% ARG bctype : type of boundary conditions
% ARG isNumeric : specifies if the solution is obtained numerically or
% analytically
% RETURN g : gradient of the parameters.
%
% SEEALSO : multiKernParamInit, multiKernCompute, sheatKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 13
    isNumeric = false;
    if nargin < 12
        bctype = 'sin';
    end
end

if strcmp(bctype, 'cos')
    cosS1 = cos(w(n+1)*s1);
    cosS2 = cos(w(m+1)*s2);
    if isNumeric
        tol = 1e-9;
        wn = (2*n+1)*pi/(2*lengthX);
        wm = (2*m+1)*pi/(2*lengthX);
        z=@(x,y)((x-y).^2).*cos(wn*x).*cos(wm*y).*exp(-((x-y).*(x-y))/(sigmax^2));
        dC = dblquad(z, 0, lengthX, 0, lengthX, tol);
        dC = (2/sigmax^3)*dC;
    else
        Wox = pwz1(n+1) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n+1)*lengthX)*pwz2(n+1);
        dWox = ((sigmax*gamma(n+1)^2)/2)*pwz1(n+1) - gamma(n+1)/(sqrt(pi)) ...
            - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n+1)/2)^2 - (lengthX/sigmax)^2))...
            *exp(-(lengthX/sigmax)^2)*exp(-gamma(n+1)*lengthX)*pwz2(n+1) + ...
            (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n+1)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            *exp(-gamma(n+1)*lengthX);        
        if n==m
            dC =  (sqrt(pi)*lengthX/2)*(real(Wox) ...
                - imag(Wox)*((sigmax^2*(2*n+1)*pi)/(4*lengthX^2))) ...
                + (sigmax*sqrt(pi)*lengthX/2)*(real(dWox) ...
                - imag(dWox)*((sigmax^2*(2*n+1)*pi)/(4*lengthX^2)) ...
                - imag(Wox)*((sigmax*(2*n+1)*pi)/(2*lengthX^2))) ...
                + sigmax*(exp(-(lengthX/sigmax)^2)*cos((2*n+1)*pi/2) - 1)  ...
                + ((lengthX^2)/sigmax)*exp(-(lengthX/sigmax)^2)*cos((2*n+1)*pi/2);
        else
            Woxm = pwz1(m+1) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m+1)*lengthX)*pwz2(m+1);
            Woxn = Wox;
            dWoxm = ((sigmax*gamma(m+1)^2)/2)*pwz1(m+1) - gamma(m+1)/(sqrt(pi)) ...
                - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(m+1)/2)^2 - (lengthX/sigmax)^2))...
                *exp(-(lengthX/sigmax)^2)*exp(-gamma(m+1)*lengthX)*pwz2(m+1) + ...
                (2/(sigmax*sqrt(pi)))*(sigmax*gamma(m+1)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
                *exp(-gamma(m+1)*lengthX);
            dWoxn = dWox;
            %             dWoxn = ((sigmax*gamma(n)^2)/2)*pwz1(n) - gamma(n)/(sqrt(pi)) ...
            %                 - ((2*lengthX^2/sigmax^3) + (2/sigmax)*((sigmax*gamma(n)/2)^2 - (lengthX/sigmax)^2))...
            %                 *exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n) + ...
            %                 (2/(sigmax*sqrt(pi)))*(sigmax*gamma(n)/2 - lengthX/sigmax)*exp(-(lengthX/sigmax)^2) ...
            %                 *exp(-gamma(n)*lengthX);
            if mod(n+m,2)==0
                dC = (lengthX/(2*sqrt(pi)*(m-n)))*(imag(Woxm) - imag(Woxn)) ...
                    + ((sigmax*lengthX)/(2*sqrt(pi)*(m-n)))*(imag(dWoxm) - imag(dWoxn));
            else
                dC = (lengthX/(2*sqrt(pi)*(m+n+1)))*(imag(Woxm) + imag(Woxn)) ...
                    + ((sigmax*lengthX)/(2*sqrt(pi)*(m+n+1)))*(imag(dWoxm) + imag(dWoxn));
            end
        end
    end
    g = dC*sum(sum((cosS1*cosS2').*covGrad));
else
    sinS1 = sin(w(n)*s1);
    sinS2 = sin(w(m)*s2);
    if isNumeric
        tol = 1e-9;
        wn = n*pi/lengthX;
        wm = m*pi/lengthX;
        z=@(x,y)((x-y).^2).*sin(wn*x).*sin(wm*y).*exp(-((x-y).*(x-y))/(sigmax^2));
        dC = dblquad(z, 0, lengthX, 0, lengthX, tol);
        dC = (2/sigmax^3)*dC;
    else
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
    end
    g = dC*sum(sum((sinS1*sinS2').*covGrad));
end
