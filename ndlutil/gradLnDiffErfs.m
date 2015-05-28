function [dlnPart, m] = gradLnDiffErfs(x1, x2, fact1, fact2),

% GRADLNDIFFERFS Compute the gradient of the log difference of two erfs.
%
%
% COPYRIGHT : Antti Honkela, 2007
  
% NDLUTIL

m = min(x1.^2, x2.^2);
dlnPart = 2/sqrt(pi) * (exp(-x1.^2 + m) .* fact1 ...
			- exp(-x2.^2 + m) .* fact2);
