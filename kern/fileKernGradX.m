function gX = fileKernGradX(kern, X, X2)

% FILEKERNGRADX Gradient of FILE kernel with respect to a point x.
% This command makes no sense for the FILE kernel. 
%
% SEEALSO fileKernParamInit, kernGradX, fileKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
