function g = biasKernGradientK(kern, x, covGrad)

% BIASKERNGRADIENTK Gradient of bias kernel wrt its parameters.

g{1} = ones(size(x, 1));
