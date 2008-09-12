function y = gpSubspaceOut(model,x)

y = NaN.*ones(size(x,1),length(model.dim));
y(:,find(model.dim)) = gpOut(model,x);

return;