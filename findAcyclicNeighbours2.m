function [S, neighbours] = findAcyclicNeighbours2(Y, k)
  
% FINDACYCLICNEIGHBOURS2 find the k nearest neighbours for each point in Y preventing cycles in the graph.
% FORMAT
% DESC returns the indices of the k nearest neighbours to each point in
% the given data matrix Y.
% ARG y : the data in which neighbours need to be found.
% ARG k : the number of neighbours that need to be found.
% RETURN ind : the indices of each points neighbours.
% RETURN D : the squared distance to each of the neighbours.
%
% COPYRIGHT : Neil D. Lawrence, 2010
%
% SEEALSO : lleOptimise, fmvuOptimise, isomapCreate

% MLTOOLS
  
  [neighboursInd, A] = findNeighbours(Y, k);
  N = size(Y, 1);
  W = spalloc(N, N, 2*size(Y, 1)*k);
  for i = 1:N
    for j = 1:k
      W(i, neighboursInd(i, j)) = -1;
      W(neighboursInd(i, j), i) = -1;
    end
  end
  L = W;
  jitter = 1e-6;
  L(1:N+1:end) = -sum(W)+jitter;
  %P = amd(L);
  %Y = Y(P, :);
  P = 1:size(L, 1);
  [UT, p, S] = chol(L, 'lower', 'vector');
  %UT = chol(L(P, P))';
  UT(1:N+1:end)=0;
  neighbours{1} = [];
  for i = 1:N-1
    neighbours{i} = find(UT(:, i));
  end
end  
