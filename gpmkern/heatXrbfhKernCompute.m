function [K, sK] = heatXrbfhKernCompute(heatKern, rbfhKern, x1, x2)

% HEATXRBFHKERNCOMPUTE Cross kernel between a HEAT and a RBF kernels.
% FORMAT
% DESC computes cross kernel terms between a HEAT kernel and a RBF kernel
% for the multiple output kernel.
% ARG heatKern : the kernel structure associated with the HEAT kernel.
% ARG rbfKern : the kernel structure associated with the RBFH kernel.
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix 
%
% FORMAT
% DESC computes cross kernel terms between a HEAT kernel and a RBF kernel
% for the multiple output kernel.
% ARG heatKern : the kernel structure associated with the HEAT kernel.
% ARG rbfKern : the kernel structure associated with the RBFH kernel.
% ARG x1 : row inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% ARG x2 : column inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% RETURN k : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix 
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
        K = zeros(length(t1)*length(s1), length(t2)*length(s2));        
    else        
        isPointwise = true;
        K = zeros(length(t1), length(t2));        
    end
else
    isPointwise = false;
    K = zeros(length(t1)*length(s1), length(t2)*length(s2));
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
   error('Not implemented yet')
else
    if isPointwise
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            Kt = simXrbfKernCompute(heatKern.sim, rbfhKern.rbf, t1, t2);
            Ks = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, i);
            K = K + Kt.*Ks;
        end
    else
        for i=1:nterms
            heatKern.sim.decay = beta(i);
            Kt = simXrbfKernCompute(heatKern.sim, rbfhKern.rbf, t1, t2);
            Ks = srbfhKernCompute(sigmax, lengthX, s1, s2, w, gamma, i);
            K = K + kron(Kt,Ks);
        end
    end
    sK = cK*K;
    K = heatKern.sensitivity*sK; 
end




