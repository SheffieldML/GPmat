function gT = simKernDiagGradX(kern, t);

% SIMKERNDIAGGRADX Gradient of SIM kernel's diagonal with respect to the
% input times t.
% FORMAT
% DESC computes the gradient of the diagonal of the single input motif
% kernel matrix with respect to the elements of the design matrix given
% in t.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG t : the input data in the form of a design matrix.
% RETURN gT : the gradients of the diagonal with respect to each element
% of t. The returned matrix has the same dimensions as t.
%
% SEEALSO : simKernParamInit, kernDiagGradX, simkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : David Luengo, 2009

% KERN

if size(t, 2) > 1
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);
t = t - kern.delay;
halfSigmaD = 0.5*sigma*kern.decay;

gT = zeros(size(t));

if (kern.isStationary == false)
    lnPart1 = lnDiffErfs(halfSigmaD, halfSigmaD-t/sigma);
    gT = kern.variance * exp(halfSigmaD*halfSigmaD - 2*kern.decay*t + lnPart1);
else
    gT = zeros(size(t));
end

if ~isfield(kern, 'isNormalised') || (kern.isNormalised == false)
    gT = gT * sigma * sqrt(pi);
end

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  error('simKerDiagGradX not implemented for gaussianInitial')
end
