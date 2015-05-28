function model = gpTimeDynamicsSetLatentValues(model, X)

% GPTIMEDYNAMICSSETLATENTVALUES Set the latent values inside the model.
%
%	Description:
%
%	MODEL = GPTIMEDYNAMICSSETLATENTVALUES(MODEL, X) Distributes GP-LVM
%	latent positions throughout the GP time dynamics model.
%	 Returns:
%	  MODEL - the updated model with the relevant latent positions in
%	   place.
%	 Arguments:
%	  MODEL - the model in which the latent positions are to be placed.
%	  X - the latent positions to be placed in the model.
%	
%	
%
%	See also
%	GPTIMEDYNAMICSCREATE, GPTIMEDYNAMICSLOGLIKELIHOOD


%	Copyright (c) 2006 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008


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

