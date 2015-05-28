function X = rocholMultiply(ch, Y);

% ROCHOLMULTIPLY Multiply by the rank one Cholesky.

% ROCHOL

%/~
% This would be the long way of doing it.
%L = rocholExtract(ch);
%X1 = L*Y;
%~/
X = zeros(size(Y));
t = zeros(1, size(Y, 2));
X(1, :) = Y(1, :)*ch.s(1);
for i = 2:ch.n
  t = t + ch.s(i-1)*ch.u(i-1)*Y(i-1, :);
  X(i, :) = ch.s(i)*Y(i, :)+ ch.v(i)*t;
end

%/~
%disp(max(max(X - X1)));
%~/

