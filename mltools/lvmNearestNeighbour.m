function err = lvmNearestNeighbour(model, lbls)

% LVMNEARESTNEIGHBOUR Give the number of errors in latent space for 1 nearest neighbour.
% FORMAT
% DESC computes the number errors for 1 nearest neighbour in latent
% space.
% ARG model : the model for which the computation is required.
% ARG lbls : the labels of the data.
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2006, 2008
  
% MLTOOLS

d = dist2(model.X, model.X);
for i = 1:size(model.X, 1); 
  d(i, i) = inf; 
end

for i= 1:size(lbls, 1); 
  lbls2(i, :) =  find(lbls(i, :));
end
[void, ind] = min(d);
err = size(model.X, 1) - sum(lbls2(ind) == lbls2);
