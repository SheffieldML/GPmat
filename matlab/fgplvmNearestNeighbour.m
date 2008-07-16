function err = fgplvmNearestNeighbour(model, lbls)

% FGPLVMNEARESTNEIGHBOUR Give the number of errors in latent space for 1 nearest neighbour.

% FGPLVM

d = dist2(model.X, model.X);
for i = 1:size(model.X, 1); 
  d(i, i) = inf; 
end

for i= 1:size(lbls, 1); 
  lbls2(i, :) =  find(lbls(i, :));
end
[void, ind] = min(d);
err = size(model.X, 1) - sum(lbls2(ind) == lbls2);
