function g = gplvmgradient(params, Y, X)

% GPLVMGRADIENT Gradient of likelihood with respect to parameters.

global SQRT % Tells us to optimise in sqrt of theta space

dataDim = size(Y, 2);
numData = size(Y, 1);

lntheta = params(end-2:end);
if nargin < 3
  latentDim = (length(params) -3)/numData;
  % We must have to optimise X as well -- it isn't given!
  X = reshape(params(1:numData*latentDim), numData, latentDim);
  xOpt = 1;
else
  xOpt = 0;
end
    
if SQRT
  theta = lntheta*lntheta;
else
  theta = exp(lntheta);
end


theta = thetaConstrain(theta);

[K, invK] = computeKernel(X, theta);

dL_dK = -dataDim/2*invK + .5*invK*Y*Y'*invK;
dk = zeros(1, 3);
[dK{1}, dK{2}] = kernelDiffParams(X, theta);
for i = 1:2
  dk(i) = sum(sum(dL_dK.*dK{i}));
end
dk(3) = -sum(diag(dL_dK)/theta(end).^2);

if SQRT
  % This line is if one uses sqrt of kernel parameters
  gParam = dk.*2.*sqrt(theta)-2*(theta.^(-1/2));
else
  % This line if you parameterise via logs
  gParam = dk.*theta-1;
end

if xOpt
  dL_dx = zeros(numData, latentDim);
  
  for q = 1:latentDim
    Kpart = kernelDiffX(X, theta, q);
    dL_dx(:, q) = 2*sum(dL_dK.*Kpart, 2) - diag(dL_dK).*diag(Kpart);
  end
  
  % look for MAP solution (weight decay)
  gX= dL_dx(:)'- X(:)';
  g = -[gX(:)' gParam];
else
  g = -gParam;
end