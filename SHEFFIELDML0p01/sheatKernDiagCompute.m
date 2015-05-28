function k = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, pwz1, pwz2, n, m)

% SHEATKERNDIAGCOMPUTE Compute a diagonal for the SHEAT kernel matrix.
%
%	Description:
%
%	K = SHEATKERNDIAGCOMPUTE(SIGMAX, LENGTHX, S, W, GAMMA, PWZ1, PWZ2,
%	N, M) computes the diagonal for a kernel matrix computed using a
%	SHEAT kernel. It gives a partial result needed for the computation
%	of the diagonal of the HEAT kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  SIGMAX - length-scale of the spatial gp prior.
%	  LENGTHX - length of the spatial domain
%	  S - inputs for which diagonal of the kernel is to be computed.
%	  W - precomputed constant.
%	  GAMMA - precomputed constant.
%	  PWZ1 - precomputed constant.
%	  PWZ2 - precomputed constant.
%	  N - integer indicating first series term
%	  M - integer indicating second series term
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, HEATKERNDIAGCOMPUTE


%	Copyright (c) 2010 Mauricio A. Alvarez


sinS1 = sin(w(n)*s);
sinS2 = sin(w(m)*s);

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

constSinS1 = const*sinS1;
k = constSinS1.*sinS2;
