function [K, const] = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, pwz1, pwz2, n, m, bctype, isNumeric)

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
% ARG bctype : type of boundary conditions
% ARG isNumeric : specifies if the solution is obtained numerically or
% analytically
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, heatKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 12
    isNumeric = false;
    if nargin < 11
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
        z=@(x,y)cos(wn*x).*cos(wm*y).*exp(-((x-y).*(x-y))/(sigmax^2));
        const = dblquad(z, 0, lengthX, 0, lengthX, tol);
    else
        Wox = pwz1(n+1) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n+1)*lengthX)*pwz2(n+1);
        if n == m            
            const = (sigmax*sqrt(pi)*lengthX/2)*(real(Wox) ...
                - imag(Wox)*((sigmax^2*(2*n+1)*pi)/(4*lengthX^2))) ...
                +(sigmax^2/2)*(exp(-(lengthX/sigmax)^2)*cos((2*n+1)*pi/2) - 1);
        else
            Woxm = pwz1(m+1) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m+1)*lengthX)*pwz2(m+1); 
            Woxn = Wox;
            if mod(n+m,2)==0               
                const = (1/(m-n))*(imag(Woxm) - imag(Woxn));                
            else
                const = (1/(m+n+1))*(imag(Woxm) + imag(Woxn));
            end
            const = (sigmax*lengthX)/(2*sqrt(pi))*const;
        end        
    end
    constCosS1 = const*cosS1;
    K = constCosS1*cosS2';    
else
    sinS1 = sin(w(n)*s1);
    sinS2 = sin(w(m)*s2);
    if isNumeric
        tol = 1e-6;
        wn = n*pi/lengthX;
        wm = m*pi/lengthX;
        z=@(x,y)sin(wn*x).*sin(wm*y).*exp(-((x-y).*(x-y))/(sigmax^2));
        const = dblquad(z, 0, lengthX, 0, lengthX, tol);
    else
        if n == m
            Wox = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
            const = (sigmax*sqrt(pi)*lengthX/2)*(real(Wox) ...
                - imag(Wox)*((sigmax^2*n*pi)/(2*lengthX^2) + (1/(n*pi)))) ...
                +(sigmax^2/2)*(exp(-(lengthX/sigmax)^2)*cos(n*pi) - 1);
        else
            if mod(n+m,2)==1
                const = 0;
            else
                Woxm = pwz1(m) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m);
                Woxn = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
                const = ((sigmax*lengthX)/(sqrt(pi)*(m^2-n^2)))*(n*imag(Woxm) - m*imag(Woxn));
            end
        end
    end    
    constSinS1 = const*sinS1;
    K = constSinS1*sinS2';    
end
