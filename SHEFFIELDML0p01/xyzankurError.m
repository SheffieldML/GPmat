function [e E] = xyzankurError(X1,X2,type)

% XYZANKURERROR Computes the error between two poses in xyz format
%
%	Description:
%
%	E = XYZANKURERROR(X1, X2) Computes the error between two poses in
%	xyz format data.
%	 Returns:
%	  E - mean error for the sequence
%	 Arguments:
%	  X1 - first set of poses
%	  X2 - second set of poses
%	
%	
%
%	See also
%	XYZANKURDRAW, XYZANKURMODIFY


%	Copyright (c) Andreas Damianou  and Neil Lawrence, 2011 Carl Henrik Ek


if(nargin<2)
    error('Too few arguments');
end
if(size(X1)~=size(X2))
    error('Dimensions mismatch');
end

E = zeros(size(X1,1),1);
for(i = 1:1:size(X1,1))
    E(i) = mean(sqrt(sum((xyzankur2joint(X1(i,:)) - xyzankur2joint(X2(i,:))).^2,2)));
end
e = mean(E);

return
