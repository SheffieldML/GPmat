function model = robTwoDynamicsCreate(q, d, latentVals)

% ROBTWODYNAMICSCREATE Create the dynamics model. 

% FGPLVM

if(q~=2)
  error('Robot Two Dynamics designed for 2-D latent spaces.');
end
if(d~=q)
  error('Input dimension must equal output dimension');
end
model.a = 20;
model.sigma2One = (pi/64).^2;
model.sigma2Two = (pi/8).^2;
model.mixThetaOne = 0.6;
model.mixThetaTwo = 0.35;
model.mixR = 0.5;

model = robTwoDynamicsSetLatentValues(model, latentVals);

% set b (the scale of the gamma distributions) so that a/b is the average jump.
aveR = mean(model.r);
model.b = model.a/aveR;
model.type = 'robTwoDynamics';

