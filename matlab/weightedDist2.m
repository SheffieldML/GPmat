% WEIGHTEDDIST2 Calculates the square of the euclidean distance between a vector and a
% matrix, where the sub-distances in every dimension are weighted by a
% given vector. 

% DESC DST = WEIGHTEDDIST2(v, X, w) Calculates the square of the euclidean
% distance between a 1xQ-dimensional vector v and a NxQ matrix X,
% with the 1xQ-dimensional vector w scaling the sub-distance
% for each of the Q dimensions. For a weight vector of ones, this function
% is equivalent to dist2.

% COPYRIGHT: Andreas C. Damianou, 2011
% SEEALSO: dist2
% SHEFFIELDML

function dst = weightedDist2(v, X, w)

N = size(X,1);
A = repmat(v, N,1);
W = repmat(w, N,1); 

%dst=sum( ((A-X).^2).*W,2)';
dst=sum( ((W.*(A-X)).^2),2)';



