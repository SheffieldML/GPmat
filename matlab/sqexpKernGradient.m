function g = sqexpKernGradient(kern, x, covGrad)

% SQEXPKERNGRADIENT Gradient of squared exponential kernel's parameters.

% IVM


[k, rbfPart, dist2xx] = sqexpKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx));
g(2) =  sum(sum(covGrad.*rbfPart))/kern.rbfVariance;
g(3) =  trace(covGrad);
g(4) =  sum(sum(covGrad));

factorVector = [kern.inverseWidth kern.rbfVariance ...
                kern.whiteVariance kern.biasVariance];
if kern.linearBound
  factors = gradFactLinearBound(factorVector);
  g = g.*factors;
else
  g = g.*factorsVector;
end
