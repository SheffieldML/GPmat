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
% COPYRIGHT : Neil D. Lawrence and Carl Henrik Ek, 2006
%
% MODIFICATIONS : Carl Henrik Ek, 2008

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

if(isfield(model,'indexIn')&&~isempty(model.indexIn))
  model.X = X(ind_in,model.indexIn);
else
  model.X = X(ind_in, :);
end
  
if(isfield(model,'indexOut')&&~isempty(model.indexOut))
  if(model.diff)
    model.y = X(ind_out,model.indexOut) - X(ind_in,model.indexOut);
  else
    model.y = X(ind_out,model.indexOut);
  end
else
  if model.diff
    model.y = X(ind_out, :) - X(ind_in, :);
  else
    model.y = X(ind_out, :);
  end
end
model.m = gpComputeM(model);

