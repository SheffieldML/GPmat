function [S, neighbours] = findAcyclicNeighbours2(Y, k)

% FINDACYCLICNEIGHBOURS2 find the k nearest neighbours for each point in Y preventing cycles in the graph.
%
%	Description:
%
%	[IND, D] = FINDACYCLICNEIGHBOURS2(Y, K) returns the indices of the k
%	nearest neighbours to each point in the given data matrix Y.
%	 Returns:
%	  IND - the indices of each points neighbours.
%	  D - the squared distance to each of the neighbours.
%	 Arguments:
%	  Y - the data in which neighbours need to be found.
%	  K - the number of neighbours that need to be found.
%	
%
%	See also
%	LLEOPTIMISE, FMVUOPTIMISE, ISOMAPCREATE


%	Copyright (c) 2010 Neil D. Lawrence

  
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
