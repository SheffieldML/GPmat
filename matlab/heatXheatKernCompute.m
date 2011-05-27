function [K, sK, sKIC] = heatXheatKernCompute(heatKern1, heatKern2, x1, x2)

% HEATXHEATKERNCOMPUTE Compute a cross kernel between two HEAT kernels.
% FORMAT
% DESC computes cross kernel terms between two HEAT kernels for
% the multiple output kernel.
% ARG heatKern1 : the kernel structure associated with the first HEAT
% kernel.
% ARG heatKern2 : the kernel structure associated with the second HEAT
% kernel.
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
% RETURN sKIC : unscaled kernel matrix associated to the initial conditions
%
% FORMAT
% DESC computes cross kernel terms between two HEAT kernels for
% the multiple output kernel.
% ARG heatKern1 : the kernel structure associated with the first HEAT
% kernel.
% ARG heatKern2 : the kernel structure associated with the second HEAT
% kernel.
% ARG x1 : row inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% ARG x2 : column inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% RETURN k : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
% RETURN sKIC : unscaled kernel matrix associated to the initial conditions
%
% SEEALSO : multiKernParamInit, multiKernCompute, heatKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    x2 = x1;
end
if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
    error('Input can only have two columns');
end
if (heatKern1.inverseWidthTime ~= heatKern2.inverseWidthTime) || ...
        (heatKern1.inverseWidthSpace ~= heatKern2.inverseWidthSpace)
    error('Kernels cannot be cross combined if they have different inverse widths.')
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
    us1 = unique(s1);
    us2 = unique(s2);
    if (length(ut1)*length(us1) == length(t1)) && ...
            (length(ut2)*length(us2) == length(t2))
        t1 = ut1; s1 = us1; t2 = ut2; s2 = us2;
        isPointwise = false;
        K = zeros(length(t1)*length(s1), length(t2)*length(s2));
        if heatKern1.includeIC
            sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
            sKIC = zeros(length(t1)*length(s1), length(t2)*length(s2));
        end
    else
        isPointwise = true;
        K = zeros(length(t1), length(t2));
        if heatKern1.includeIC
            sK = zeros(length(t1), length(t2));
            sKIC = zeros(length(t1), length(t2));
        end
    end
else
    isPointwise = false;
    K = zeros(length(t1)*length(s1), length(t2)*length(s2));
    if heatKern1.includeIC
        sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
        sKIC = zeros(length(t1)*length(s1), length(t2)*length(s2));
    end
end

% Although this is done in heatKernExpandParam.m, we do it here again as a
% precaution.

heatKern1.sim.inverseWidth = heatKern1.inverseWidthTime;
heatKern2.sim.inverseWidth = heatKern1.inverseWidthTime;

sigmax = sqrt(2/heatKern1.inverseWidthSpace);
lengthX = heatKern1.lengthX;
nterms = heatKern1.nTerms;
decay1 = heatKern1.decay;
diff1 = heatKern1.diffusion;
decay2 = heatKern2.decay;
diff2 = heatKern2.diffusion;


if strcmp(heatKern1.pde, 'cos')
    w = ((2*(0:(nterms-1))+1)*(pi/(2*lengthX)))';
    gamma = sqrt(-1)*w;
    beta1 = decay1 + diff1*(w.^2);
    beta2 = decay2 + diff2*(w.^2);
    z1 = sigmax*gamma/2;
    z2 = lengthX/sigmax + z1;
    wz1 = wofzPoppe(sqrt(-1)*z1);
    wz2 = wofzPoppe(sqrt(-1)*z2);
    cK = 4/(lengthX^2);
    if isPointwise
        for i=0:nterms-1
            for j=0:nterms-1                
                heatKern1.sim.decay = beta1(i+1);
                heatKern2.sim.decay = beta2(j+1);
                Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, heatKern1.pde);
                K = K + Kt.*Ks;                
            end
        end
    else
        for i=0:nterms-1
            for j=0:nterms-1                
                heatKern1.sim.decay = beta1(i+1);
                heatKern2.sim.decay = beta2(j+1);
                Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, heatKern1.pde);
                K = K + kron(Kt,Ks);                
            end
        end
    end
    sK = cK*K;
    K = heatKern1.sensitivity*heatKern2.sensitivity*sK;
    if nargout >2
        sKIC = 0;
    end    
else    
    % Precompute some terms
    w = ((1:nterms)*(pi/lengthX))';
    gamma = sqrt(-1)*w;
    beta1 = decay1 + diff1*(w.^2);
    beta2 = decay2 + diff2*(w.^2);
    z1 = sigmax*gamma/2;
    z2 = lengthX/sigmax + z1;
    wz1 = wofzPoppe(sqrt(-1)*z1);
    wz2 = wofzPoppe(sqrt(-1)*z2);
    cK = 4/(lengthX^2);
    
    if heatKern1.includeIC
        if (heatKern1.inverseWidthSpaceIC ~= heatKern2.inverseWidthSpaceIC)
            error('Kernels cannot be cross combined if they have different inverse widths for the initial conditions.')
        end
        sigmah = sqrt(2/heatKern1.inverseWidthSpaceIC);
        z1h = sigmah*gamma/2;
        z2h = lengthX/sigmah + z1h;
        wz1h = wofzPoppe(sqrt(-1)*z1h);
        wz2h = wofzPoppe(sqrt(-1)*z2h);
        if isPointwise
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        KtIC = exp(-beta1(i)*t1)*(exp(-beta2(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        sKIC = sKIC + KtIC.*KsIC;
                        sK = sK + Kt.*Ks;
                    end
                end
            end
        else
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        KtIC = exp(-beta1(i)*t1)*(exp(-beta2(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        sKIC = sKIC + kron(KtIC,KsIC);
                        sK = sK + kron(Kt,Ks);
                    end
                end
            end
        end
        sKIC = cK*sKIC;
        sK = cK*sK;
        K = (heatKern1.sensitivity*heatKern2.sensitivity)*sK ...
            + (heatKern1.sensitivityIC*heatKern2.sensitivityIC)*sKIC;
    else
        if isPointwise
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        K = K + Kt.*Ks;
                    end
                end
            end
        else
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        K = K + kron(Kt,Ks);
                    end
                end
            end
        end
        sK = cK*K;
        K = heatKern1.sensitivity*heatKern2.sensitivity*sK;
        if nargout >2
            sKIC = 0;
        end
    end    
end


K = real(K);
sK = real(sK);




