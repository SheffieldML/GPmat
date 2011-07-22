function [g1, g2] = heatXrbfhKernGradient(heatKern, rbfhKern, x1, x2, covGrad)

% HEATXRBFHKERNGRADIENT Gradient wrt parameters between a HEAT and a RBFH.
% FORMAT
% DESC computes the gradients wrt parameters of a cross kernel term between 
% a HEAT kernel and a RBFH kernel for the multiple output kernel.
% ARG heatKern : the kernel structure associated with the HEAT kernel.
% ARG rbfKern : the kernel structure associated with the RBFH kernel.
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
% DESC computes the gradients wrt parameters of a cross kernel term between 
% a HEAT kernel and a RBFH kernel for the multiple output kernel.
% ARG heatKern : the kernel structure associated with the HEAT kernel.
% ARG rbfKern : the kernel structure associated with the RBFH kernel.
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
if (heatKern.inverseWidthTime ~= rbfhKern.inverseWidthTime) || ...
        (heatKern.inverseWidthSpace ~= rbfhKern.inverseWidthSpace)
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
    else
        isPointwise = true;
        sK = zeros(length(t1), length(t2));      
    end
else
    isPointwise = false;
    sK = zeros(length(t1)*length(s1), length(t2)*length(s2));
end
% Although this is done in heatKernExpandParam.m, we do it here again as a
% precaution.

heatKern.sim.inverseWidth = heatKern.inverseWidthTime;
rbfhKern.rbf.inverseWidth = rbfhKern.inverseWidthTime;

sigmax = sqrt(2/heatKern.inverseWidthSpace);
lengthX = heatKern.lengthX;
nterms = heatKern.nTerms;
decay = heatKern.decay;
diff = heatKern.diffusion;

% Precompute some terms
w = ((1:nterms)*(pi/lengthX))';
gamma = sqrt(-1)*w;
beta = decay + diff*(w.^2);
cK = 2/lengthX;

if heatKern.includeIC
   error('Not implemented yet.')    
else
    g1 = zeros(1,5);
    g2 = zeros(1,2);
    if isPointwise
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            Kt = simXrbfKernCompute(heatKern.sim, rbfhKern.rbf, t1, t2);
            [Ks, wz1, wz2] = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, i);
            covGradt = covGrad.*Ks;
            covGrads = covGrad.*Kt;
            g1A = simXrbfKernGradient(heatKern.sim, rbfhKern.rbf, t1, t2, covGradt);
            g1(1) = g1(1) + g1A(1);          % Decay for heat kernel
            g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the heat kernel
            g1(3) = g1(3) + g1A(2);          % Inverse width for time
            g1B = srbfhKernGradient(sigmax, lengthX, s1, s2, w, gamma, i, covGrads, wz1, wz2);
            g1B = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
            g1(4) = g1(4) + g1B;
            sK = sK + Kt.*Ks;
        end
    else
        covGradt = zeros(length(t1), length(t2));
        covGrads = zeros(length(s1), length(s2));
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            Kt = simXrbfKernCompute(heatKern.sim, rbfhKern.rbf, t1, t2);
            [Ks, wz1, wz2] = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, i);
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
            g1A = simXrbfKernGradient(heatKern.sim, rbfhKern.rbf, t1, t2, covGradt);
            g1(1) = g1(1) + g1A(1);          % Decay for heat kernel
            g1(2) = g1(2) + (w(i)^2)*g1A(1); % Diffusion rate for the heat kernel
            g1(3) = g1(3) + g1A(2);          % Inverse width for time
            g1B = srbfhKernGradient(sigmax, lengthX, s1, s2, w, gamma, i, covGrads, wz1, wz2);
            g1B = -(1/sqrt(2*heatKern.inverseWidthSpace^3))*g1B; % Transforms to the derivative of the inverse width
            g1(4) = g1(4) + g1B;
            sK = sK + kron(Kt, Ks);
        end
    end
    g1(1:4) = heatKern.sensitivity*cK*g1(1:4);
    g1(5) = cK*(sum(sum(covGrad.*sK))); 
end




