function k = disimKernDiagCompute(kern, t)

% DISIMKERNDIAGCOMPUTE Compute diagonal of DISIM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the driven input
%  single input motif kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : disimKernParamInit, kernDiagCompute, kernCreate, disimKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2009

% KERN

if size(t, 2) > 1 
  error('Input can only have one column');
end

l = sqrt(2/kern.inverseWidth);
t = t;
delta = kern.di_decay;
D = kern.decay;
halfLD = 0.5*l*D;
halfLDelta = 0.5*l*delta;

[lnPart1, signs1] = lnDiffErfs(halfLDelta - t/l, ...
			       halfLDelta);
[lnPart2, signs2] = lnDiffErfs(halfLDelta + t/l, ...
			       halfLDelta);

lnCommon = halfLDelta .^ 2 -(D+delta)*t - log(2*delta) - log(D-delta);
lnFact2 = (D+delta)*t - log(D + delta);


if abs(D-delta) < .1,
  h = signs1 .* exp(lnCommon + lnPart1) ...
      .* ((exp((D-delta)*t) - 1) / (D - delta) + 1/(D+delta)) ...
      + signs2 .* exp(lnCommon + lnFact2 + lnPart2);
else
  lnFact1a = (D - delta) * t + log(D + delta) - log(D^2 - delta^2);
  lnFact1b = log(2*delta) - log(D^2 - delta^2);

  h = signs1 .* exp(lnCommon + lnFact1a + lnPart1) ...
      - signs1 .* exp(lnCommon + lnFact1b + lnPart1) ...
      + signs2 .* exp(lnCommon + lnFact2 + lnPart2);
end

[lnPart1p, signs1p] = lnDiffErfs(halfLD - t/l, ...
				 halfLD);
[lnPart2p, signs2p] = lnDiffErfs(halfLD + t/l, ...
				 halfLD);

lnCommonp = halfLD.^2 - 2*D*t - log(delta^2 - D^2);
lnFact2p = 2*D*t - log(2*D);

if abs(D-delta) < .1,
  hp = signs1p .* exp(lnCommonp + lnPart1p) ...
       .* ((exp((D-delta)*t) - 1) / (D - delta) + 1/(2*D)) ...
       + signs2p .* exp(lnCommonp + lnFact2p + lnPart2p);
else
  lnFact1ap = log(D + delta) - log(delta - D) - log(2*D);
  lnFact1bp = (D-delta)*t - log(delta - D);

  hp = signs1p .* exp(lnCommonp + lnFact1ap + lnPart1p) ...
       - signs1p .* exp(lnCommonp + lnFact1bp + lnPart1p) ...
       + signs2p .* exp(lnCommonp + lnFact2p + lnPart2p);
end

k = 2*real(h+hp);
k = 0.5*k*sqrt(pi)*l;
k = kern.rbf_variance*kern.di_variance*kern.variance*k;

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  k = k + kern.initialVariance * kern.variance * ...
      ((exp(-kern.di_decay*t) - exp(-kern.decay*t)) ./ (kern.decay-kern.di_decay)).^2;
end
