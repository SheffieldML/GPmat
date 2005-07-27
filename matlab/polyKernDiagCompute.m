function k = polyKernDiagCompute(kern, x)

% POLYKERNDIAGCOMPUTE Compute diagonal of polynomial kernel.

% KERN

k =  kern.variance*(sum(x.*x, 2)*kern.weightVariance + kern.biasVariance).^kern.degree;
