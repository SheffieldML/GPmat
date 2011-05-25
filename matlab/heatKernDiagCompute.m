function [k, sk, skIC] = heatKernDiagCompute(heatKern, x)

% HEATKERNDIAGCOMPUTE Diagonal of a kernel matrix for a HEAT kernel.
% DESC computes the diagonal of the kernel matrix for the HEAT kernel
% given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
% RETURN sk : unscaled version of the diagonal.
% RETURN sk : unscaled version of the diagonal for the kernel matrix of the
% initial conditions.
%
% SEEALSO : heatKernParamInit, kernDiagCompute, kernCreate, heatKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(x, 2) ~= 2
    error('Input can only have two columns');
end

option = 'cos';

if strcmp(option, 'cos')    
    [k, sk] = heatXheatKernCompute(heatKern, heatKern, x, x);    
    k = diag(k);
    sk = diag(sk);
else
    % Split the domain into time domain and spatial domain and account for
    % missing values. If there are no missing values the computation of the
    % kernel is a pointwise prodruct, otherwise it is a kronecker product.
    t = x(x(:,1)~=Inf,1);
    s = x(x(:,2)~=Inf,2);
    
    if (length(t) == length(s))
        ut = unique(t);
        us = unique(s);
        if (length(ut)*length(us) == length(t))
            t = ut;
            s = us;
            isPointwise = false;
            sk = zeros(length(t)*length(s), 1);
            if heatKern.includeIC
                skIC = zeros(length(t)*length(s), 1);
            end
        else
            isPointwise = true;
            sk = zeros(length(t), 1);
            if heatKern.includeIC
                skIC = zeros(length(t), 1);
            end
        end
    else
        isPointwise = false;
        sk = zeros(length(t)*length(s), 1);
        if heatKern.includeIC
            skIC = zeros(length(t)*length(s), 1);
        end
    end
    
    % Although this is done in heatKernExpandParam.m, we do it here again as a
    % precaution.
    
    heatKern.sim.inverseWidth = heatKern.inverseWidthTime;
    
    sigmax  = sqrt(2/heatKern.inverseWidthSpace);
    lengthX = heatKern.lengthX;
    nterms  = heatKern.nTerms;
    decay   = heatKern.decay;
    diff    = heatKern.diffusion;
    
    % Precompute some terms
    w = ((1:nterms)*(pi/lengthX))';
    gamma = sqrt(-1)*w;
    beta = decay + diff*(w.^2);
    z1 = sigmax*gamma/2;
    z2 = lengthX/sigmax + z1;
    wz1 = wofzPoppe(sqrt(-1)*z1);
    wz2 = wofzPoppe(sqrt(-1)*z2);
    cK = 4/(lengthX^2);
    
    simLocal = heatKern.sim;
    
    if heatKern.includeIC
        sigmah = sqrt(2/heatKern.inverseWidthSpaceIC);
        z1h = sigmah*gamma/2;
        z2h = lengthX/sigmah + z1h;
        wz1h = wofzPoppe(sqrt(-1)*z1h);
        wz2h = wofzPoppe(sqrt(-1)*z2h);
        if isPointwise
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                kt = simKernDiagCompute(heatKern.sim, t);
                ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
                ktIC = exp(-2*beta(i)*t);
                ksIC = sheatKernDiagCompute(sigmah, lengthX, s, w, gamma, wz1h, wz2h, i, i);
                skIC = skIC + ktIC.*ksIC;
                sk = sk + kt.*ks;
                for j=1:i-1
                    if (mod(i+j,2)==0)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        %kt = diag(simXsimKernCompute(heatKern.sim, simLocal, t, t));
                        kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                        ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                        ktIC = exp(-(beta(i)+beta(j))*t);
                        ksIC = sheatKernDiagCompute(sigmah, lengthX, s, w, gamma, wz1h, wz2h, i, j);
                        skIC = skIC + 2*ktIC.*ksIC;
                        sk = sk + 2*kt.*ks;
                    end
                end
            end
        else
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                kt = simKernDiagCompute(heatKern.sim, t);
                ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
                ktIC = exp(-2*beta(i)*t);
                ksIC = sheatKernDiagCompute(sigmah, lengthX, s, w, gamma, wz1h, wz2h, i, i);
                skIC = skIC + kron(ktIC, ksIC);
                sk = sk + kron(kt, ks);
                for j=1:i-1
                    if (mod(i+j,2)==0)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        %kt = diag(simXsimKernCompute(heatKern.sim, simLocal, t, t));
                        kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                        ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                        ktIC = exp(-(beta(i)+beta(j))*t);
                        ksIC = sheatKernDiagCompute(sigmah, lengthX, s, w, gamma, wz1h, wz2h, i, j);
                        skIC = skIC + 2*kron(ktIC, ksIC);
                        sk = sk + 2*kron(kt, ks);
                    end
                end
            end
        end
        sk = cK*sk;
        skIC = cK*skIC;
        k = (heatKern.sensitivity^2)*sk + (heatKern.sensitivityIC^2)*skIC;
    else
        if isPointwise
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                kt = simKernDiagCompute(heatKern.sim, t);
                ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
                sk = sk + kt.*ks;
                for j=1:i-1
                    if (mod(i+j,2)==0)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        %kt = diag(simXsimKernCompute(heatKern.sim, simLocal, t, t));
                        kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                        ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                        sk = sk + 2*kt.*ks;
                    end
                end
            end
        else
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                kt = simKernDiagCompute(heatKern.sim, t);
                ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
                sk = sk + kron(kt, ks);
                for j=1:i-1
                    if (mod(i+j,2)==0)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        %kt = diag(simXsimKernCompute(heatKern.sim, simLocal, t, t));
                        kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                        ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                        sk = sk + kron(2*kt, ks);
                    end
                end
            end
        end
        sk = cK*sk;
        k = (heatKern.sensitivity^2)*sk;
        if nargin > 2
            skIC = 0;
        end
    end        
end

k = real(k);
sk = real(sk);
