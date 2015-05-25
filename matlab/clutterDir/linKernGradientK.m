function G = linKernGradientK(kern, x)

% LINKERNGRADIENTK Gradient of lin kernel wrt its parameters.

linPart = linKernCompute(kern, x);
G{1} = linPart/kern.variance;
