function X = rocholBackSub(ch, Y);

% ROCHOLBACKSUB Backsubstitute the representation of the rank one Cholesky.

% ROCHOL

%/~
% This would be the long way of doing it.
%L = rocholExtract(ch);
%X1 = L'\Y;
%~/
X = zeros(size(Y));
t = zeros(1, size(Y, 2));
X(ch.n, :) = Y(ch.n, :)/ch.s(ch.n);
for i = ch.n-1:-1:1
  t = t + ch.v(i+1)*X(i+1, :);
  X(i, :) = Y(i, :)/ch.s(i) - ch.u(i)*t;
end

%/~
%disp(max(max(X - X1)));
%~/

