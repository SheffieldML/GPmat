function neighboursInd = findNeighbours(Y, k)
  
% FINDNEIGHBOURS find the k nearest neighbours for each point in Y.
% FORMAT
% DESC returns the indices of the k nearest neighbours to each point in
% the given data matrix Y.
% ARG y : the data in which neighbours need to be found.
% ARG k : the number of neighbours that need to be found.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : lleCreate, mvuCreate, isomapCreate

% MLTOOLS
  
Y2 = sum(Y.*Y, 2);
D = repmat(Y2', size(Y, 1), 1) + repmat(Y2, 1, size(Y, 1)) - 2*Y*Y';
D(1:size(D, 1)+1:end) = inf;
[void, ind] = sort(D);
neighboursInd = ind(1:k, :)';
