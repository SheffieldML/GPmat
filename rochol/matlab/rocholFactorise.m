function ch = rocholFactorise(v);

% ROCHOLFACTORISE Rank one Cholesky factorise.

% ROCHOL

% Return structure which represents factorisation of I+vv^T

if size(v, 2) ~= 1
  error('v should be a row vector');
end
t = 1;
ch.n = size(v, 1);
ch.v = v;
ch.s = zeros(size(v));
ch.u = zeros(size(v));
for i = 1:ch.n
  tnew = t + ch.v(i)*ch.v(i);
  ch.s(i) = sqrt(tnew/t);
  ch.u(i) = ch.v(i)/tnew;
  t = tnew;
end
