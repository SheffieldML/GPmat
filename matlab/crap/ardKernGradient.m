function g = ardKernGradient(kern, x, covGrad)

% ARDKERNGRADIENT Gradient of ard kernel's parameters.

% IVM

[k, rbfPart, linearPart, dist2xx] = ardKernCompute(x, kern);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx))*kern.inverseWidth;
g(2) =  sum(sum(covGrad.*rbfPart));
g(3) =  sum(sum(covGrad.*eye(size(x, 1))))*kern.whiteVariance;
g(4) =  sum(sum(covGrad))*kern.biasVariance;
g(5) =  sum(sum(covGrad.*linearPart));
for i = 1:size(x, 2)
  g(5+i) =  sum(sum(covGrad.*((kern.linearVariance)*x(:, i)*x(:, i)' ...
                                      -.5*(kern.inverseWidth)*dist2(x(:, i), ...
                                                    x(:, i)) ...
                                      .*rbfPart)))*kern.inputScales(i);
end