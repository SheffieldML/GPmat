function G = whiteKernGradientK(kern, x)

% WHITEKERNGRADIENTK Gradient of white noise kernel wrt its parameters.

G{1} = eye(size(x, 1));
