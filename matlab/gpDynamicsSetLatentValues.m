function model = gpDynamicsSetLatentValues(model, X)

% GPDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% FORMAT
% DESC Distributes GP-LVM latent positions throughout the GP dynamics model. 
% ARG model : the model in which the latent positions are to be
% placed.
% ARG X : the latent positions to be placed in the model.
% RETURN model : the updated model with the relevant latent
% positions in place.
%
% SEEALSO : gpDynamicsCreate, gpDynamicsLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence and Cark Henrik Ek, 2006

% FGPLVM

ind_in = [];
ind_out = []; 
startVal=1;
for i = 1:length(model.seq)
  endVal = model.seq(i);
  ind_in = [ind_in startVal:endVal-1];
  ind_out = [ind_out startVal+1:endVal];
  startVal = endVal + 1;
end

model.X = X(ind_in, :);
if model.diff
  model.y = X(ind_out, :) - X(ind_in, :);
else
  model.y = X(ind_out, :);
end

for i = 1:model.d
  model.m(:, i) = (model.y(:, i) - model.bias(i));
  if model.scale(i)
    model.m(:, i) = model.y(:, i)/model.scale(i);
  end
end

