function gX = fileKernGradX(kern, X, X2)

% FILEKERNGRADX Gradient of file stored kernel with respect to X (zeros).

% KERN

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
