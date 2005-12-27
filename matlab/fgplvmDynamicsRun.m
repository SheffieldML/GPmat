function fgplvmDynamicsRun

% FGPLVMDYNAMICSRUN Visualise the manifold.

% FGPLVM

global visualiseInfo

while visualiseInfo.clicked & visualiseInfo.runDynamics
  % This should be changed to a model specific visualisation.
  x = visualiseInfo.latentPos(1);
  y = visualiseInfo.latentPos(2);
  set(visualiseInfo.latentHandle, 'xdata', x, 'ydata', y);
  fhandle = str2func([visualiseInfo.model.type 'PosteriorMeanVar']);
  [mu, varsigma] = fhandle(visualiseInfo.model, visualiseInfo.latentPos);
  if isfield(visualiseInfo.model, 'noise')
    Y = noiseOut(visualiseInfo.model.noise, mu, varsigma);
  else
    Y = mu;
  end
  visualiseInfo.visualiseModify(visualiseInfo.visHandle, ...
                                Y, visualiseInfo.varargin{:});
  fhandle = str2func([visualiseInfo.model.type ...
                      'DynamicsPosteriorMeanVar']);
  [mu, var] = fhandle(visualiseInfo.model, visualiseInfo.latentPos);
  visualiseInfo.latentPos=gsamp(mu, diag(var), 1);
  pause(0.05)
end



