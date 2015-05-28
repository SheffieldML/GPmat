function D1 = dist1(X, Y)

% DIST1	Calculates absolute distance (i.e. L1 norm) between two sets of
% points. The idea of the function is similar to Ian Nabney's Netlab DIST2
% function, although the computation is not so efficient in this case.
%
% FORMAT
% DESC takes two matrices of vectors and calculates the L1 distance between
% each of the points in those matrices. Both matrices must be of the same
% column dimension.  If X has M rows and N columns, and Y has L rows and
% N columns, then the result has M rows and L columns. The (I, J)-th entry of
% D1 is the L1 distance from the I-th row of X to the j-th row of Y.
% ARG X: First N-dimensional Input vector in the form of a design matrix.
% ARG Y: Second N-dimensional Input vector in the form of a design matrix.
% RETURN D1: Pointwise L1 distance matrix between X and Y.
%
% SEEALSO: dist2 (Netlab)
%
% COPYRIGHT : David Luengo, 2009

% NDLUTIL

[Nx, dimX] = size(X);
[Ny, dimY] = size(Y);
if dimX ~= dimY
	error('Data dimension of X and Y does not match')
end

if Nx <= Ny
    for i = 1:Nx
        D1(i, :) = sum(abs(repmat(X(i,:), Ny, 1) - Y), 2)';
    end
else
    for i = 1:Ny
        D1(:, i) = sum(abs(X - repmat(Y(i,:), Nx, 1)), 2);
    end
end
