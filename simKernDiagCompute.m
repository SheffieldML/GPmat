function k = simKernDiagCompute(kern, t)

% SIMKERNDIAGCOMPUTE Compute diagonal of SIM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the single input motif kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : simKernParamInit, kernDiagCompute, kernCreate, simKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Antti Honkela, 2008
%
% MODIFICATIONS : David Luengo, 2009
%
% MODIFICATIONS : Mauricio Alvarez, 2009

% KERN

if size(t, 2) > 1 
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);
t = t - kern.delay;
halfSigmaD = 0.5*sigma*kern.decay;

if (kern.isStationary == false)
    lnPart1 = lnDiffErfs(halfSigmaD + t/sigma, halfSigmaD);
    lnPart2 = lnDiffErfs(halfSigmaD, halfSigmaD - t/sigma);
    h = exp(halfSigmaD*halfSigmaD + lnPart1)...
        - exp(halfSigmaD*halfSigmaD-(2*kern.decay*t) + lnPart2);
else
    lnPart1 = lnDiffErfs(inf, halfSigmaD);
    h = exp(halfSigmaD*halfSigmaD + lnPart1) * ones(size(t));
end

if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    k = (kern.sensitivity*kern.sensitivity)*h/(2*kern.decay);
else
    k = kern.variance*h/(2*kern.decay);
end

if ~isfield(kern, 'isNormalised') || (kern.isNormalised == false)
    k = sqrt(pi)*sigma*k;
end

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  k = k + kern.initialVariance*exp(-2*kern.decay*t);
end
