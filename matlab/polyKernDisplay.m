function polyKernDisplay(kern)

% POLYKERNDISPLAY Display parameters of polynomial kernel.

% KERN

fprintf('Polynomial kernel Variance: %2.4f\n', kern.variance)
fprintf('Polynomial weight variance: %2.4f\n', kern.weightVariance)
fprintf('Polynomial bias variance: %2.4f\n', kern.biasVariance)