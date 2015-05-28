function X = rocholTransMultiply(ch, Y);

% ROCHOLTRANSMULTIPLY Multiply by the transposed version of the rank one Cholesky.
% ROCHOL

%/~
% This would be the long way of doing it.
%L = rocholExtract(ch);
%X1 = L'*Y;
%~/
if size(Y, 1) ~= ch.n
  error('Inner matrix dimensions do not match');
end
X = zeros(size(Y));
t = zeros(1, size(Y, 2));
X(ch.n, :) = Y(ch.n, :)*ch.s(ch.n);
for i = ch.n-1:-1:1
  t = t + ch.v(i+1)*Y(i+1, :);
  X(i, :) = ch.s(i)*(Y(i, :)+ ch.u(i)*t);
end

%/~
%disp(max(max(X - X1)));
%~/

