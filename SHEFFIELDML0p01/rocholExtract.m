function L = rocholExtract(ch)

% ROCHOLEXTRACT Extract the lower triangular matrix from the Cholesky structure.
%
%	Description:
%	L = rocholExtract(ch)
%

L = diag(ch.s);
for j = 1:ch.n
  for i = j+1:ch.n
    L(i, j) = ch.v(i)*ch.s(j)*ch.u(j);
  end
end
