 function model = oldIvmRemovePoint(model, i)

% IVMREMOVEPOINT Remove the least informative point from the IVM.


index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set']);
end

model.beta(i) = 0;
model.m(i) = 0;

vector = zeros(length(model.I), 1);
vector(index) = 1e20;
[L, D, Ltilde, Dprime, Lprime, beta, p] = updateCholesky(model.Sigma.L, [], vector, '+');

LtildeInv = eye(size(Ltilde))/Ltilde;
%model.Sigma.M = diag(1./sqrt(diag(Dprime)))*LtildeInv*sqrt(D)*model.Sigma.M;

model.Sigma.M = diag(1./sqrt(diag(Dprime)))*forwardSub(sqrt(D)*model.Sigma.M, ...
						 beta, p);


model.Sigma.M(index, :) = [];
model.kern.Kstore(:, index) = [];

model.Sigma.L = Lprime*diag(sqrt(diag(Dprime)));
model.Sigma.L(index, :) = [];
model.Sigma.L(:, index) =[];

model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;
model.I(index) = [];
model.J = [model.J i];
model.J = sort(model.J);

model.varSigma = model.kern.diagK - sum(model.Sigma.M.*model.Sigma.M)';
Linv =eye(size(model.Sigma.L))/model.Sigma.L;

model.mu =  model.kern.Kstore*Linv'*Linv*model.m(model.I);

model = feval([model.noise.type 'UpdateParams'], model, model.J);