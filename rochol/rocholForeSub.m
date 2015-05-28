function X = rocholForeSub(ch, Y);

% ROCHOLFORESUB Foreward substitute the representation of the rank one Cholesky.

% ROCHOL

%/~
% This would be the long way of doing it.
%L = rocholExtract(ch);
%X2 = L\Y;
%~/

X = zeros(size(Y));
t = zeros(1, size(Y, 2));
X(1, :) = Y(1, :)/ch.s(1);
for i = 2:ch.n
  if i == 2
    t = Y(1, :)*ch.u(1);
  else
    t = t + ch.s(i-1)*ch.u(i-1)*X(i-1, :);
  end
  X(i, :) = (Y(i, :)-ch.v(i)*t)/ch.s(i);
end

%/~
%disp(max(max(X - X2)));
%~/
