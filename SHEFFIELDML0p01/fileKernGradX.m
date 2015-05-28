function gX = fileKernGradX(kern, X, X2)

% FILEKERNGRADX Gradient of FILE kernel with respect to a point x.
%
%	Description:
%	This command makes no sense for the FILE kernel.
%	
%	
%
%	See also
%	% SEEALSO FILEKERNPARAMINIT, KERNGRADX, FILEKERNDIAGGRADX


%	Copyright (c) 2005, 2006 Neil D. Lawrence



gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
