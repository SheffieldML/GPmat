function model = downdateM(model, ni)

% DOWNDATEM 


index = find(model.I == ni);
if isempty(index)
  error(['Point ' num2str(ni) ' is not in active set']);
end

sigmavec = model.kern.Kstore(:, index) - model.Sigma.M'*model.Sigma.M(:, ni);
nuprime = 1/(1/model.beta(ni) - model.varSigma(ni));
%/~
oldVarSigma = model.varSigma;
if any(model.varSigma<0)
  warning('Variance less than zero')
end
%~/
model.varSigma = model.varSigma + nuprime*sigmavec.*sigmavec;
%/~
if any(model.varSigma<0)
  warning('Variance less than zero')
end
%~/
gprime = nuprime*(model.mu(ni)-model.m(ni));
model.mu = model.mu + gprime*sigmavec;

model.m(ni) = 0;
model.beta(ni) = 0;

model.kern.Kstore(:, index) = [];

model.I(index) = [];
model.Sigma.L = chol(model.kern.Kstore(model.I, :) + diag(model.beta(model.I)))';
model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;
model.Sigma.M = model.Sigma.Linv*model.kern.Kstore';
model.J = [model.J ni];
model.J = sort(model.J);

return

v = model.Sigma.Linv(index:end, index);
b = zeros(size(v));
d = zeros(size(v));
t(1) = v(1)*v(1);
for i = 2:length(v)
  t(i) = t(i-1)+ v(i)*v(i);
  d(i) = sqrt(t(i)/t(i-1));
  b(i) = v(i)/t(i);
end
b = b(2:end);
d = d(2:end);
v = v(2:end);

% downdate model.Sigma.L
% downdate model.
for i = 1:length(v)
  Ltilde(i, i) = d(i);
  lambda(i) = d(i);
  for j = 1:i-1
    a = v(i)*b(j)*d(j);
    Ltilde(i, j) = a;
    lambda(i) = lambda(i) + a;
    
  end
end

% m = model.Sigma.M(index, :);
% if ~isempty(b)
%   model.Sigma.M(index+1:end, :) = forwardSub(model.Sigma.M(index+1:end, :), v, b, d);
% end
% for i = 1:length(v)
%   model.Sigma.M(index+i, :) = model.Sigma.M(index+i, :) - m*lambda(i);
% end
% model.Sigma.M(index, :) = [];
LtildeFull= zeros(size(model.Sigma.L));
for i = 1:index
  LtildeFull(i, i) = 1;
end
LtildeFull(index:end, index) = inf;
LtildeFull(index+1:end, index+1:end) = Ltilde;

[L, D, LtildeFull2, Dprime, Lprime, beta, p] = updateCholesky(model.Sigma.L, [], 1e10*model.Sigma.Linv(:, index), '+');

% l = model.Sigma.Linv(index, :);
% if ~isempty(b)
%   model.Sigma.Linv(index+1:end, :) = forwardSub(model.Sigma.Linv(index+1:end, :), v, b, d);
% end
% for i = 1:length(v)
%   model.Sigma.Linv(index+i, :) = model.Sigma.Linv(index+i, :) - l*lambda(i);
% end
% model.Sigma.Linv(index, :) = [];
% model.Sigma.Linv(:, index) = [];

% model.Sigma.L = Lprime*sqrt(Dprime);
% model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;
% model.Sigma.L(index, :) = [];
% model.Sigma.L(:, index) = [];
% model.Sigma.Linv(index, :) = [];
% model.Sigma.Linv(:, index) = [];
% model.kern.Kstore(:, index) = [];
% model.Sigma.M = model.Sigma.Linv*model.kern.Kstore';
% L = model.Sigma.L(index+1:end, index+1:end);
% for i =1:length(v)
%   rho = L(i, i)*v(i); 
%   model.Sigma.L(i+index, i+index) = L(i, i);
%   sigma = 0;
%   for j = i-1:-1:1
%     sigma = sigma + rho;
%     rho = L(i, j)*v(j);
%     model.Sigma.L(i+index, j) = L(i, j) + sigma*b(j);
%   end
% end
% for i = 1:length(v)
%   model.Sigma.L(i+index, i+index:end) = model.Sigma.L(i+index, i+index:end)*d(i);
% end
% model.Sigma.L(index, :) = [];
% model.Sigma.L(:, index) = [];
%model.Sigma.L = eye(size(model.Sigma.Linv))\model.Sigma.Linv;

