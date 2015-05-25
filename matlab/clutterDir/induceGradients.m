function g = induceGradients(xVals, model)

% INDUCEGRADIENTS Gradients of the log det A wrt inducing variables

X_u = reshape(xVals, model.k, model.q);
K_uu = kernCompute(model.kern, X_u);
K_uf = kernCompute(model.kern, X_u, model.X);
gA = pdinv(1/model.beta*K_uu + K_uf*K_uf');
gK_u = gA*1/model.beta;
gK_uf = 2*gA*K_uf;

%%% Compute Gradients with respect to X_u %%%
gKX = kernGradX(model.kern, X_u, X_u);

% The 2 accounts for the fact that covGrad is symmetric
gKX = gKX*2;
dgKX = kernDiagGradX(model.kern, X_u);
for i = 1:model.k
  gKX(i, :, i) = dgKX(i, :);
end

% Allocate space for gX_u
gX_u = zeros(model.k, model.q);
% Compute portion associated with gK_u
for i = 1:model.k
  for j = 1:model.q
    gX_u(i, j) = gKX(:, j, i)'*gK_u(:, i);
  end
end


% Compute portion associated with gK_uf
gKX_uf = kernGradX(model.kern, X_u, model.X);
for i = 1:model.k
  for j = 1:model.q
    gX_u(i, j) = gX_u(i, j) + gKX_uf(:, j, i)'*gK_uf(i, :)';
  end
end
g = -gX_u(:)';


