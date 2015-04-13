function n2 = ndlDist2(X1, X2)
% DIST2	Calculates squared distance between two sets of points.
% FORMAT
% DESC this rewrite of the NETLAB dist2 takes two matrices of vectors and calculates the
%	squared Euclidean distance between them.  Both matrices must be of
%	the same column dimension.  If X has M rows and N columns, and C has
%	L rows and N columns, then the result has M rows and L columns.  The
%	I, Jth entry is the  squared distance from the Ith row of X to the
%	Jth row of C.
% ARG X1 : first matrix for distances.
% ARG X2 : second matrix for distances.
% RETURN n2 : the distances between the matrices.
%
% SEEALSO : dist2
%
% COPYRIGHT : Neil D. Lawrence, 2008

% NDLUTIL
  
  
if nargin <2
  numData = size(X1, 1);
  d2 = 2*repmat(sum(X1.*X1, 2)', numData, 1) ...
       - 2*X1*X1';
else
  
  numData1 = size(X1, 1);
  numData2 = size(X2, 2);
  
  d2 = repmat(sum(X1.*X1, 2)', numData2, 1) ...
       - 2*X1*X2' ...
       + repmat(sum(X2.*X2, 2)', numData1, 1);
end