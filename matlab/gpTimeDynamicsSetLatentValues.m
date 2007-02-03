function model = gpTimeDynamicsSetLatentValues(model, X)

% GPTIMEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
% FORMAT
% DESC Distributes GP-LVM latent positions throughout the GP time dynamics model. 
% ARG model : the model in which the latent positions are to be
% placed.
% ARG X : the latent positions to be placed in the model.
% RETURN model : the updated model with the relevant latent
% positions in place.
%
% SEEALSO : gpTimeDynamicsCreate, gpTimeDynamicsLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006

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

if model.diff
  model.y = X(ind_out, :) - X(ind_in, :);
else
  model.y = X(ind_out, :);
end

model.m = gpComputeM(model);

