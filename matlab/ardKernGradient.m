function g = ardKernGradient(kern, x, covGrad)

% ARDKERNGRADIENT Gradient of ard kernel's parameters.

% IVM

[k, rbfPart, linearPart, dist2xx] = ardKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx));
g(2) =  sum(sum(covGrad.*rbfPart))/kern.rbfVariance;
g(3) =  sum(sum(covGrad.*eye(size(x, 1))));
g(4) =  sum(sum(covGrad));
g(5) =  sum(sum(covGrad.*linearPart))/kern.linearVariance;
for i = 1:size(x, 2)
  g(5+i) =  sum(sum(covGrad.*((kern.linearVariance)*x(:, i)*x(:, i)' ...
                                      -.5*(kern.inverseWidth)*dist2(x(:, i), ...
                                                    x(:, i)) ...
                                      .*rbfPart)))*kern.inputScales(i)*(1-kern.inputScales(i));
end
factorVector = [kern.inverseWidth kern.rbfVariance ...
                kern.whiteVariance kern.biasVariance kern.linearVariance];
if kern.linearBound
  factors = gradFactLinearBound(factorVector);
  g(1:5) = g(1:5).*factors;
else
  g(1:5) = g(1:5).*factorVector;
end