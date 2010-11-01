function [K, sK, sKIC] = heatKernCompute(heatKern, x1, x2)

% HEATKERNCOMPUTE Compute a kernel matrix for a HEAT kernel.
% FORMAT computes the kernel matrix for the heat kernel function given
% inputs associated with rows and columns.
% ARG heatKern : the kernel structure associated with the HEAT
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
% RETURN sKIC : unscaled kernel matrix associated to the initial conditions
%
% FORMAT
% DESC computes the kernel matrix for the single input motif
% kernel given a design matrix of inputs.
% ARG heatKern : the kernel structure associated with HEATvkernel.
% RETURN k : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
% RETURN sKIC : unscaled kernel matrix associated to the initial conditions
%
% SEEALSO : multiKernParamInit, multiKernCompute, heatKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 3
    x2 = x1;
    if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
        error('Input can only have two columns');
    end
    % Split the domain into time domain and spatial domain and account for
    % missing values. If there are no missing values the computation of the
    % kernel is a pointwise prodruct, otherwise it is a kronecker product.
    t1 = x1(x1(:,1)~=Inf,1);
    t2 = x2(x2(:,1)~=Inf,1);
    s1 = x1(x1(:,2)~=Inf,2);
    s2 = x2(x2(:,2)~=Inf,2);

    
    
    
    if (length(t1) == length(s1)) && (length(t2) == length(s2))
        ut1 = unique(t1);
        ut2 = unique(t2);
        if (length(ut1) == length(t1)) || (length(ut2) == length(t2))
            isPointwise = true;
            K = zeros(length(t1), length(t2));
            if heatKern.includeIC
                sK = zeros(length(t1), length(t2));
                sKIC = zeros(length(t1), length(t2));
            end
        else
            us1 = unique(s1);us2 = unique(s2);
            t1 = ut1; s1 = us1; t2 = ut2; s2 = us2;
            isPointwise = false;
            K = zeros(length(t1)*length(s1), length(t2)*length(s2));
            if heatKern.includeIC
                sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
                sKIC = zeros(length(t1)*length(s1), length(t2)*length(s2));
            end            
        end
    else
        isPointwise = false;
        K = zeros(length(t1)*length(s1), length(t2)*length(s2));
        if heatKern.includeIC
            sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
            sKIC = zeros(length(t1)*length(s1), length(t2)*length(s2));
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
                Kt = simKernCompute(heatKern.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, i);
                KtIC = exp(-beta(i)*t1)*(exp(-beta(i)*t2)');
                KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, i);
                sKIC = sKIC + KtIC.*KsIC;
                sK = sK + Kt.*Ks;
                for j=1:i-1
                    if (mod(i+j,2)==1)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        Kt = simXsimKernCompute(heatKern.sim, simLocal, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        KtIC = exp(-beta(i)*t1)*(exp(-beta(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        tempo = KtIC.*KsIC;
                        sKIC = sKIC + tempo + tempo';
                        tempo = Kt.*Ks;
                        sK = sK + tempo + tempo';
                    end
                end
            end
        else
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                Kt = simKernCompute(heatKern.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, i);
                KtIC = exp(-beta(i)*t1)*(exp(-beta(i)*t2)');
                KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, i);
                sKIC = sKIC + kron(KtIC,KsIC);
                sK = sK + kron(Kt,Ks);
                for j=1:i-1
                    if (mod(i+j,2)==1)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        Kt = simXsimKernCompute(heatKern.sim, simLocal, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        KtIC = exp(-beta(i)*t1)*(exp(-beta(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        tempo = kron(KtIC, KsIC);
                        sKIC = sKIC + tempo + tempo';
                        tempo = kron(Kt, Ks);
                        sK = sK + tempo + tempo';
                    end
                end
            end
        end
        sK = cK*sK;
        sKIC = cK*sKIC;
        K = (heatKern.sensitivity^2)*sK + (heatKern.sensitivityIC^2)*sKIC;
    else
        if isPointwise
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                Kt = simKernCompute(heatKern.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, i);
                K = K + Kt.*Ks;
                for j=1:i-1
                    if (i == j) || (mod(i+j,2)==1)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        Kt = simXsimKernCompute(heatKern.sim, simLocal, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        tempo = Kt.*Ks;
                        K = K + tempo + tempo';
                    end
                end
            end
        else
            for i=1:nterms
                heatKern.sim.decay = beta(i);
                Kt = simKernCompute(heatKern.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, i);
                K = K + kron(Kt, Ks);
                for j=1:i-1
                    if (i == j) || (mod(i+j,2)==1)
                        heatKern.sim.decay = beta(i);
                        simLocal.decay = beta(j);
                        Kt = simXsimKernCompute(heatKern.sim, simLocal, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        tempo = kron(Kt, Ks);
                        K = K + tempo + tempo';
                    end
                end
            end
        end
        sK = cK*K;
        K = (heatKern.sensitivity^2)*sK;
        if nargout >2
            sKIC = 0;
        end
    end
else
    [K, sK, sKIC] = heatXheatKernCompute(heatKern, heatKern, x1, x2);
end





