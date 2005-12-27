function model = fgplvmAddDynamics(model, kernType, snRatio)

% FGPLVMADDDYNAMICS Add a dynamics kernel to the model.

% FGPLVM
% Set up dynamics kernel.
model.dynamics.q = model.q;
model.dynamics.d = model.q;
model.dynamics.N = model.N -1;
model.dynamics.approx = model.approx;
model.dynamics.k = model.k;

if ~strcmp(model.approx, 'ftc')
  model.dynamics.sigma2 = 1e-6;
  model.dynamics.sigma2Transform = 'negLogLogit';
end

if nargin< 3 
  snRatio = 100;
end
if isstruct(kernType)
  model.dynamics.kern = kernType;
else
  model.dynamics.kern = kernCreate(model.X(1:end-1, :), kernType);
  signalSize = kernGetVariance(model.dynamics.kern);
  for i = 1:length(model.dynamics.kern.comp)
    model.dynamics.kern.comp{i}.variance = 0.01/signalSize;
  end
  
  model.dynamics.kern = kernSetWhite(model.dynamics.kern, signalSize/(snRatio*snRatio));
end

% If the approximation is pitc, need to add block ends to dynamics.
switch model.approx
 case 'pitc'
  model.dynamics.blockEnd = model.blockEnd;
  model.dynamics.blockEnd(end) = model.dynamics.blockEnd(end)-1;
  if model.dynamics.blockEnd(end) == model.dynamics.blockEnd(end-1)
    model.dynamics.blockEnd(end) = [];
  end
end
params = fgplvmExtractParam(model);
model = fgplvmExpandParam(model, params);
