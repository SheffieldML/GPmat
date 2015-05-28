function [g1, g2] = heatXheatKernGradient(heatKern1, heatKern2, x1, x2, covGrad)

% HEATXHEATKERNGRADIENT Gradient wrt parameters between two HEAT kernels.
% FORMAT
% DESC computes the gardients wrt parameters of a cross kernel term between
% two HEAT kernels for the multiple output kernel.
% ARG heatKern1 : the kernel structure associated with the first HEAT
% kernel.
% ARG heatKern2 : the kernel structure associated with the second HEAT
% kernel.
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
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
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
%
% SEEALSO : heatXheatKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 5
    covGrad = x2;
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
        sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
        if heatKern1.includeIC
            sKIC = zeros(length(t1)*length(s1), length(t2)*length(s2));
        end
    else
        isPointwise = true;
        sK = zeros(length(t1), length(t2));
        if heatKern1.includeIC
            sKIC = zeros(length(t1), length(t2));
        end
    end
else
    isPointwise = false;
    sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
    if heatKern1.includeIC
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
    g1 = zeros(1,5);
    g2 = zeros(1,5);
    if isPointwise
        for i=0:nterms-1
            for j=0:nterms-1                
                heatKern1.sim.decay = beta1(i+1);
                heatKern2.sim.decay = beta2(j+1);
                Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, heatKern1.pde);
                covGradt = covGrad.*Ks;
                covGrads = covGrad.*Kt;
                [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                g1(1) = g1(1) + g1A(1);   % Decay for first kernel
                g1(2) = g1(2) + (w(i+1)^2)*g1A(1); % Diffusion rate for the first kernel
                g1(3) = g1(3) + g1A(2);          % Inverse width for time
                g2(1) = g2(1) + g2A(1);          % Decay for second kernel
                g2(2) = g2(2) + (w(j+1)^2)*g2A(1); % Diffusion rate for the second kernel
                g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads, heatKern1.pde);
                g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                g1(4) = g1(4) + g1B;
                sK = sK + Kt.*Ks;               
            end
        end
    else
        covGradt = zeros(length(t1), length(t2));
        covGrads = zeros(length(s1), length(s2));
        for i=0:nterms-1
            for j=0:nterms-1               
                heatKern1.sim.decay = beta1(i+1);
                heatKern2.sim.decay = beta2(j+1);
                Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, heatKern1.pde);
                % These loops might slow everything
                startOne = 1;
                endOne = 0;
                for k=1:length(t1)
                    endOne = endOne + length(s1);
                    startTwo = 1;
                    endTwo = 0;
                    for l=1:length(t2)
                        endTwo = endTwo + length(s2);
                        covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
                        startTwo = endTwo + 1;
                    end
                    startOne = endOne + 1;
                end
                for k=1:length(s1)
                    indRows = (k:length(s1):(length(t1)*length(s1)))';
                    for l=1:length(s2)
                        indCols = l:length(s2):(length(t2)*length(s2));
                        covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
                    end
                end
                [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                g1(1) = g1(1) + g1A(1);   % Decay for first kernel
                g1(2) = g1(2) + (w(i+1)^2)*g1A(1); % Diffusion rate for the first kernel
                g1(3) = g1(3) + g1A(2);          % Inverse width for time
                g2(1) = g2(1) + g2A(1);          % Decay for second kernel
                g2(2) = g2(2) + (w(j+1)^2)*g2A(1); % Diffusion rate for the second kernel
                g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads, heatKern1.pde);
                g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                g1(4) = g1(4) + g1B;
                sK = sK + kron(Kt, Ks);                
            end
        end
    end
    g1(1:4) = heatKern1.sensitivity*heatKern2.sensitivity*cK*g1(1:4);
    g2(1:4) = heatKern1.sensitivity*heatKern2.sensitivity*cK*g2(1:4);
    g1(5) = heatKern2.sensitivity*cK*(sum(sum(covGrad.*sK)));
    g2(5) = heatKern1.sensitivity*cK*(sum(sum(covGrad.*sK)));
    
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
        g1 = zeros(1,7);
        g2 = zeros(1,7);
        sigmah = sqrt(2/heatKern1.inverseWidthSpaceIC);
        z1h = sigmah*gamma/2;
        z2h = lengthX/sigmah + z1h;
        wz1h = wofzPoppe(sqrt(-1)*z1h);
        wz2h = wofzPoppe(sqrt(-1)*z2h);
        cKf = heatKern1.sensitivity*heatKern2.sensitivity*cK;
        cKh = heatKern1.sensitivityIC*heatKern2.sensitivityIC*cK;
        if isPointwise
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        covGradt = covGrad.*Ks;
                        covGrads = covGrad.*Kt;
                        [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                        g1(1) = g1(1) + cKf*g1A(1);   % Decay for first kernel
                        g1(2) = g1(2) + cKf*(w(i)^2)*g1A(1); % Diffusion rate for the first kernel
                        g1(3) = g1(3) + cKf*g1A(2);          % Inverse width for time
                        g2(1) = g2(1) + cKf*g2A(1);          % Decay for second kernel
                        g2(2) = g2(2) + cKf*(w(j)^2)*g2A(1); % Diffusion rate for the second kernel
                        g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads);
                        g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                        g1(4) = g1(4) + cKf*g1B;
                        KtIC = exp(-beta1(i)*t1)*(exp(-beta2(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        covGradtIC = covGrad.*KsIC;
                        covGradsIC = covGrad.*KtIC;
                        g1AIC = sum(sum(((- t1.*exp(-beta1(i)*t1))*(exp(-beta2(j)*t2)')).*covGradtIC));
                        g2AIC = sum(sum((exp(-beta1(i)*t1)*((-t2.*exp(-beta2(j)*t2))')).*covGradtIC));
                        g1(1) = g1(1) + cKh*g1AIC;    % Decay for first kernel
                        g1(2) = g1(2) + cKh*(w(i)^2)*g1AIC;  % Diffusion rate for the first kernel
                        g2(1) = g2(1) + cKh*g2AIC;           % Decay for second kernel
                        g2(2) = g2(2) + cKh*(w(j)^2)*g2AIC;  % Diffusion rate for the second kernel
                        g1BIC = sheatKernGradient(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j, covGradsIC);
                        g1BIC = -(1/sqrt(2*heatKern1.inverseWidthSpaceIC^3))*g1BIC; % Transforms to the derivative of the inverse width
                        g1(6) = g1(6) + cKh*g1BIC;
                        sKIC = sKIC + KtIC.*KsIC;
                        sK = sK + Kt.*Ks;
                    end
                end
            end
        else
            covGradt = zeros(length(t1), length(t2));
            covGrads = zeros(length(s1), length(s2));
            covGradtIC = zeros(length(t1), length(t2));
            covGradsIC = zeros(length(s1), length(s2));
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        % These loops might slow everything
                        startOne = 1;
                        endOne = 0;
                        for k=1:length(t1)
                            endOne = endOne + length(s1);
                            startTwo = 1;
                            endTwo = 0;
                            for l=1:length(t2)
                                endTwo = endTwo + length(s2);
                                covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
                                startTwo = endTwo + 1;
                            end
                            startOne = endOne + 1;
                        end
                        for k=1:length(s1)
                            indRows = (k:length(s1):(length(t1)*length(s1)))';
                            for l=1:length(s2)
                                indCols = l:length(s2):(length(t2)*length(s2));
                                covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
                            end
                        end
                        [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                        g1(1) = g1(1) + cKf*g1A(1);   % Decay for first kernel
                        g1(2) = g1(2) + cKf*(w(i)^2)*g1A(1); % Diffusion rate for the first kernel
                        g1(3) = g1(3) + cKf*g1A(2);          % Inverse width for time
                        g2(1) = g2(1) + cKf*g2A(1);          % Decay for second kernel
                        g2(2) = g2(2) + cKf*(w(j)^2)*g2A(1); % Diffusion rate for the second kernel
                        g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads);
                        g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                        g1(4) = g1(4) + cKf*g1B;
                        KtIC = exp(-beta1(i)*t1)*(exp(-beta2(j)*t2)');
                        KsIC = sheatKernCompute(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j);
                        % These loops might slow everything
                        startOne = 1;
                        endOne = 0;
                        for k=1:length(t1)
                            endOne = endOne + length(s1);
                            startTwo = 1;
                            endTwo = 0;
                            for l=1:length(t2)
                                endTwo = endTwo + length(s2);
                                covGradtIC(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*KsIC));
                                startTwo = endTwo + 1;
                            end
                            startOne = endOne + 1;
                        end
                        for k=1:length(s1)
                            indRows = (k:length(s1):(length(t1)*length(s1)))';
                            for l=1:length(s2)
                                indCols = l:length(s2):(length(t2)*length(s2));
                                covGradsIC(k,l) = sum(sum(covGrad(indRows, indCols).*KtIC));
                            end
                        end
                        g1AIC = sum(sum(((- t1.*exp(-beta1(i)*t1))*(exp(-beta2(j)*t2)')).*covGradtIC));
                        g2AIC = sum(sum((exp(-beta1(i)*t1)*((-t2.*exp(-beta2(j)*t2))')).*covGradtIC));
                        g1(1) = g1(1) + cKh*g1AIC;    % Decay for first kernel
                        g1(2) = g1(2) + cKh*(w(i)^2)*g1AIC;  % Diffusion rate for the first kernel
                        g2(1) = g2(1) + cKh*g2AIC;           % Decay for second kernel
                        g2(2) = g2(2) + cKh*(w(j)^2)*g2AIC;  % Diffusion rate for the second kernel
                        g1BIC = sheatKernGradient(sigmah, lengthX, s1, s2, w, gamma, wz1h, wz2h, i, j, covGradsIC);
                        g1BIC = -(1/sqrt(2*heatKern1.inverseWidthSpaceIC^3))*g1BIC; % Transforms to the derivative of the inverse width
                        g1(6) = g1(6) + cKh*g1BIC;
                        sKIC = sKIC + kron(KtIC, KsIC);
                        sK = sK + kron(Kt, Ks);
                    end
                end
            end
        end
        g1(5) = heatKern2.sensitivity*cK*(sum(sum(covGrad.*sK)));
        g2(5) = heatKern1.sensitivity*cK*(sum(sum(covGrad.*sK)));
        g1(7) = heatKern2.sensitivityIC*cK*(sum(sum(covGrad.*sKIC)));
        g2(7) = heatKern1.sensitivityIC*cK*(sum(sum(covGrad.*sKIC)));
    else
        g1 = zeros(1,5);
        g2 = zeros(1,5);
        if isPointwise
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        covGradt = covGrad.*Ks;
                        covGrads = covGrad.*Kt;
                        [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                        g1(1) = g1(1) + g1A(1);   % Decay for first kernel
                        g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the first kernel
                        g1(3) = g1(3) + g1A(2);          % Inverse width for time
                        g2(1) = g2(1) + g2A(1);          % Decay for second kernel
                        g2(2) = g2(2) + (w(j)^2)*g2A(1); % Diffusion rate for the second kernel
                        g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads);
                        g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                        g1(4) = g1(4) + g1B;
                        sK = sK + Kt.*Ks;
                    end
                end
            end
        else
            covGradt = zeros(length(t1), length(t2));
            covGrads = zeros(length(s1), length(s2));
            for i=1:nterms
                for j=1:nterms
                    if (i == j) || (mod(i+j,2)==0)
                        heatKern1.sim.decay = beta1(i);
                        heatKern2.sim.decay = beta2(j);
                        Kt = simXsimKernCompute(heatKern1.sim, heatKern2.sim, t1, t2);
                        Ks = sheatKernCompute(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j);
                        % These loops might slow everything
                        startOne = 1;
                        endOne = 0;
                        for k=1:length(t1)
                            endOne = endOne + length(s1);
                            startTwo = 1;
                            endTwo = 0;
                            for l=1:length(t2)
                                endTwo = endTwo + length(s2);
                                covGradt(k,l) = sum(sum(covGrad(startOne:endOne, startTwo:endTwo).*Ks));
                                startTwo = endTwo + 1;
                            end
                            startOne = endOne + 1;
                        end
                        for k=1:length(s1)
                            indRows = (k:length(s1):(length(t1)*length(s1)))';
                            for l=1:length(s2)
                                indCols = l:length(s2):(length(t2)*length(s2));
                                covGrads(k,l) = sum(sum(covGrad(indRows, indCols).*Kt));
                            end
                        end
                        [g1A, g2A] = simXsimKernGradient(heatKern1.sim, heatKern2.sim, t1, t2, covGradt);
                        g1(1) = g1(1) + g1A(1);   % Decay for first kernel
                        g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the first kernel
                        g1(3) = g1(3) + g1A(2);          % Inverse width for time
                        g2(1) = g2(1) + g2A(1);          % Decay for second kernel
                        g2(2) = g2(2) + (w(j)^2)*g2A(1); % Diffusion rate for the second kernel
                        g1B = sheatKernGradient(sigmax, lengthX, s1, s2, w, gamma, wz1, wz2, i, j, covGrads);
                        g1B = -(1/sqrt(2*heatKern1.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
                        g1(4) = g1(4) + g1B;
                        sK = sK + kron(Kt, Ks);
                    end
                end
            end
        end
        g1(1:4) = heatKern1.sensitivity*heatKern2.sensitivity*cK*g1(1:4);
        g2(1:4) = heatKern1.sensitivity*heatKern2.sensitivity*cK*g2(1:4);
        g1(5) = heatKern2.sensitivity*cK*(sum(sum(covGrad.*sK)));
        g2(5) = heatKern1.sensitivity*cK*(sum(sum(covGrad.*sK)));
    end
end
g1 = real(g1);
g2 = real(g2);


