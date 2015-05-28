function g = heatKernDiagGradient(heatKern, x, covDiag)

% HEATKERNDIAGGRADIENT Gradient of the HEAT kernel's diagonal wrt parameters.
% FORMAT
% DESC computes the gradient of functions of the diagonal of the
% heat kernel matrix with respect to the parameters of the kernel. The
% parameters' gradients are returned in the order given by the
% heatKernExtractParam command.
% ARG kern : the kernel structure for which the gradients are
% computed.
% ARG x : the input times for which the gradient is being computed.
% ARG factors : partial derivatives of the function of interest with
% respect to the diagonal elements of the kernel.
% RETURN g : gradients of the relevant function with respect to each
% of the parameters. Ordering should match the ordering given in
% heatKernExtractParam.
%
% SEEALSO : heatKernParamInit, kernDiagGradient, heatKernExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if size(x, 2) ~= 2
    error('Input can only have two columns');
end

% Split the domain into time domain and spatial domain and account for
% missing values. If there are no missing values the computation of the
% kernel is a pointwise prodruct, otherwise it is a kronecker product.
t = x(x(:,1)~=Inf,1);
s = x(x(:,2)~=Inf,2);

if (length(t) == length(s))
    ut = unique(t);
    us = unique(s);
    if (length(ut)*length(us) == length(t))
        s = us;
        t = ut;
        isPointwise = false;
        sk = zeros(length(t)*length(s), 1);                
    else
        isPointwise = true;        
        sk = zeros(length(t), 1);
    end
else
    isPointwise = false;
    sk = zeros(length(t)*length(s), 1);
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
    error('Not implemented yet')
else
    g = zeros(1,5);
    if isPointwise
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            kt = simKernDiagCompute(heatKern.sim, t);
            ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
            covDiagt = covDiag.*ks;
            covDiags = covDiag.*kt;
            gA = simKernDiagGradient(heatKern.sim, t, covDiagt);
            g(1) = g(1) + gA(1);            % Decay for first kernel
            g(2) = g(2) + (w(i)^2)*gA(1);   % Diffusion rate for the first kernel
            g(3) = g(3) + gA(2);            % Inverse width for time
            gB = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i, covDiags);
            gB = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*gB; % Transforms to the derivative of the inverse width
            g(4) = g(4) + gB;
            sk = sk + kt.*ks;
            for j=1:i-1
                if (mod(i+j,2)==0)
                    heatKern.sim.decay = beta(i);
                    simLocal.decay = beta(j);
                    kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                    ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                    covDiagt = covDiag.*ks;
                    covDiags = covDiag.*kt;
                    [g1A, g2A] = simXsimKernDiagGradient(heatKern.sim, simLocal, t, covDiagt);
                    g(1) = g(1) + 2*(g1A(1) + g2A(1));                      % Decay
                    g(2) = g(2) + 2*((w(i)^2)*g1A(1) + (w(j)^2)*g2A(1));    % Diffusion
                    g(3) = g(3) + 2*g1A(2);                               % Inverse width for time
                    gB = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j, covDiags);
                    gB = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*gB; % Transforms to the derivative of the inverse width
                    g(4) = g(4) + 2*gB;
                    sk = sk + 2*kt.*ks;
                end
            end
        end
    else
        covDiagt = zeros(length(t), 1);
        covDiags = zeros(length(s), 1);
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            kt = simKernDiagCompute(heatKern.sim, t);
            ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i);
            startOne = 1;
            endOne = 0;
            for k=1:length(t)
                endOne = endOne + length(s);
                covDiagt(k) = sum(covDiag(startOne:endOne).*ks);
                startOne = endOne + 1;
            end
            for k=1:length(s)
                indRows = (k:length(s):(length(t)*length(s)))';
                covDiags(k) = sum(covDiag(indRows).*kt);
            end
            gA = simKernDiagGradient(heatKern.sim, t, covDiagt);
            g(1) = g(1) + gA(1);            % Decay for first kernel
            g(2) = g(2) + (w(i)^2)*gA(1);   % Diffusion rate for the first kernel
            g(3) = g(3) + gA(2);            % Inverse width for time
            gB = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, wz1, wz2, i, i, covDiags);
            gB = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*gB; % Transforms to the derivative of the inverse width
            g(4) = g(4) + gB;
            sk = sk + kron(kt, ks);
            for j=1:i-1
                if (mod(i+j,2)==0)
                    heatKern.sim.decay = beta(i);
                    simLocal.decay = beta(j);
                    kt = simXsimKernDiagCompute(heatKern.sim, simLocal, t);
                    ks = sheatKernDiagCompute(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j);
                    startOne = 1;
                    endOne = 0;
                    for k=1:length(t)
                        endOne = endOne + length(s);
                        covDiagt(k) = sum(covDiag(startOne:endOne).*ks);
                        startOne = endOne + 1;
                    end
                    for k=1:length(s)
                        indRows = (k:length(s):(length(t)*length(s)))';
                        covDiags(k) = sum(covDiag(indRows).*kt);
                    end
                    [g1A, g2A] = simXsimKernDiagGradient(heatKern.sim, simLocal, t, covDiagt);
                    g(1) = g(1) + 2*(g1A(1) + g2A(1));                      % Decay
                    g(2) = g(2) + 2*((w(i)^2)*g1A(1) + (w(j)^2)*g2A(1));    % Diffusion
                    g(3) = g(3) + 2*g1A(2);                               % Inverse width for time
                    gB = sheatKernDiagGradient(sigmax, lengthX, s, w, gamma, wz1, wz2, i, j, covDiags);
                    gB = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*gB; % Transforms to the derivative of the inverse width
                    g(4) = g(4) + 2*gB;
                    sk = sk + kron(2*kt, ks);
                end
            end
        end
    end
    g(1:4) = (heatKern.sensitivity^2)*cK*g(1:4);
    g(5) = 2*heatKern.sensitivity*cK*(sum(covDiag.*sk));
end




