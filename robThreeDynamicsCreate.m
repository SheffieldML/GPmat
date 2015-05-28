function model = robThreeDynamicsCreate(q, d, latentVals)

% ROBTHREEDYNAMICSCREATE Create the dynamics model. 

% FGPLVM

if(q~=2)
  error('Robot Three Dynamics designed for 2-D latent spaces.');
end
if(d~=q)
  error('Input dimension must equal output dimension');
end
model.lambda = 1;
model.sigma2 = 1/64;
model = robThreeDynamicsSetLatentValues(model, latentVals);

model.type = 'robThreeDynamics';

