function ll = lvmScoreModel(model)
  
% LVMSCOREMODEL Score model with a GP log likelihood.

% MLTOOLS
  
  options = gpOptions('ftc');
  gpmod = gpCreate(model.q, model.d, model.X, model.Y, options);
  
  gpmod = gpOptimise(gpmod);
  
  ll = modelLogLikelihood(gpmod);
  
end
