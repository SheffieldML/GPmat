function f = induceObjective(xVals, model)

% INDUCEGRADIENTS Gradients of the log det A.


X_u = reshape(xVals, model.k, model.q);
K_uu = kernCompute(model.kern, X_u);
K_uf = kernCompute(model.kern, X_u, model.X);
f = -logdet(1/model.beta*K_uu + K_uf*K_uf');
