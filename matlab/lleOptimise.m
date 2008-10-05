function model = lleOptimise(model, display, iters)

% LLEOPTIMISE Optimise an LLE model.
% FORMAT
% DESC optimises an mixtures of Gaussians model via the 
% expectation maximisation algorithm.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : mmpcaCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

model.indices = findNeighbours(model.Y, model.k);
model.W = spalloc(model.N, model.N, model.N*model.k);
for i = 1:model.N
  Ytemp = model.Y(model.indices(i, :), :);
  %  C = cov([model.Y(i, :); Ytemp]');
  Ytemp = Ytemp - repmat(model.Y(i, :), model.k, 1);
  C = Ytemp*Ytemp';
  if model.d<model.k
    %C = C + trace(C)*1e-3*eye(model.k+1);
    C = C + trace(C)*1e-3*eye(model.k);
  end
  %[fullU, v] = eig(C);
  %[void, ind] = min(diag(v));
  %what = fullU(:, ind);
  [U, jitter] = jitChol(C);
  y = U'\ones(model.k, 1);
  what = U\y;
  
  %what = (C\ones(model.k, 1));
  what = what/sum(what);
  model.W(i, model.indices(i, :)) = what';
  model.W(i, i) = -1;
  %    model.W(i, [i model.indices(i, :)]) = what';
end

options.disp = 0; 
options.isreal = 1; 
options.issym = 1; 
%[m, v] = svds(model.W', model.q+1, 0);

%model.P = speye(model.N) - model.W - model.W' +model.W'*model.W;
model.P = model.W'*model.W;
if isoctave
  warning('No eigs function in Octave');
  % Nasty hack for eigenvalue problem in Octave.
  [m, v] = eig(model.P);
  [v, order] = sort(diag(v));
  v = diag(v(1:model.q+1));
  m = m(:, order);
  m = m(:, 1:model.q+1);
else
  [m, v] = eigs_r11(model.P, model.q+1, 'sm', options);
end
[void, ind] = min(diag(v));
model.X = m(:, [1:(ind-1) (ind+1):end]);