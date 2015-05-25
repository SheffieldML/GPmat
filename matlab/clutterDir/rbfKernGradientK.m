function G = rbfKernGradientK(kern, x, covGrad)

% RBFKERNGRADIENTK Gradient of rbf kernel's parameters.

[k, dist2xx] = rbfKernCompute(kern, x);

G{1} = - .5*k.*dist2xx;
G{2} =  k/kern.variance;
