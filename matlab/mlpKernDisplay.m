function mlpKernDisplay(kern)

% MLPKERNDISPLAY Display parameters of multi-layer perceptron kernel.

% KERN

% KERN


fprintf('MLP kernel Variance: %2.4f\n', kern.variance)
fprintf('MLP weight variance: %2.4f\n', kern.weightVariance)
fprintf('MLP bias variance: %2.4f\n', kern.biasVariance)