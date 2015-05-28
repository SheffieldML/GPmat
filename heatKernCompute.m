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

% if nargin < 3
%     x2 = x1;
% end
%
% [K, sK, sKIC] = heatXheatKernCompute(heatKern, heatKern, x1, x2);
%
% K = real(K);
% sK = real(sK);
% sKIC = real(sKIC);

if nargin < 3 && ~strcmp(heatKern.pde, 'cos')
    if size(x1, 2) ~= 2 
        error('Input can only have two columns');
    end
    % Split the domain into time domain and spatial domain and account for
    % missing values. If there are no missing values the computation of the
    % kernel is a pointwise prodruct, otherwise it is a kronecker product.
    t1 = x1(x1(:,1)~=Inf,1);
    s1 = x1(x1(:,2)~=Inf,2);
    if (length(t1) == length(s1))
        ut1 = unique(t1);        
        us1 = unique(s1);        
        if (length(ut1)*length(us1) == length(t1))                
            t1 = ut1; s1 = us1; 
            isPointwise = false;
            K = zeros(length(t1)*length(s1));
        else
            isPointwise = true;
            K = zeros(length(t1));
        end
    else
        isPointwise = false;
        K = zeros(length(t1)*length(s1));
    end
    
    % Although this is done in heatKernExpandParam.m, we do it here again as a
    % precaution.
    
    heatKern.sim.inverseWidth = heatKern.inverseWidthTime;
    heatKernTempo = heatKern;    
    
    sigmax = sqrt(2/heatKern.inverseWidthSpace);
    lengthX = heatKern.lengthX;
    nterms = heatKern.nTerms;
    decay = heatKern.decay;
    diff = heatKern.diffusion;
    % Precompute some terms
    w = ((1:nterms)*(pi/lengthX))';
    gamma = sqrt(-1)*w;
    beta = decay + diff*(w.^2);    
    z1 = sigmax*gamma/2;
    z2 = lengthX/sigmax + z1;
    wz1 = wofzPoppe(sqrt(-1)*z1);
    wz2 = wofzPoppe(sqrt(-1)*z2);
    cK = 4/(lengthX^2);    
    if isPointwise
        for i=1:nterms
            heatKern.sim.decay = beta(i);            
            Kt = simXsimKernCompute(heatKern.sim, heatKern.sim, t1, t1);
            Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i);
            K = K + Kt.*Ks;            
            for j=1:i-1
                if  mod(i+j,2)==0
                    heatKern.sim.decay = beta(i);
                    heatKernTempo.sim.decay = beta(j);
                    Kt = simXsimKernCompute(heatKern.sim, heatKernTempo.sim, t1, t1);
                    Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, j);
                    K = K + Kt.*Ks + (Kt').*(Ks');
                end
            end
        end        
    else
        for i=1:nterms
            heatKern.sim.decay = beta(i);            
            Kt = simXsimKernCompute(heatKern.sim, heatKern.sim, t1, t1);
            Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, i);
            K = K + kron(Kt,Ks);            
            for j=1:i-1
                if  mod(i+j,2)==0
                    heatKern.sim.decay = beta(i);
                    heatKernTempo.sim.decay = beta(j);
                    Kt = simXsimKernCompute(heatKern.sim, heatKernTempo.sim, t1, t1);
                    Ks = sheatKernCompute(sigmax, lengthX, s1, s1, w, gamma, wz1, wz2, i, j);
                    K = K + kron(Kt,Ks) + kron(Kt', Ks');
                end
            end
        end
    end
    sK = cK*K;
    K = heatKern.sensitivity^2*sK;
    if nargout >2
        sKIC = 0;
    end    
    K = real(K);
    sK = real(sK);
else
    [K, sK, sKIC] = heatXheatKernCompute(heatKern, heatKern, x1, x2);
end



