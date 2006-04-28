function g = gpReversibleDynamicsLogLikeGradients(model)

% GPREVERSIBLEDYNAMICSLOGLIKEGRADIENTS Gradients of the GP reversible dynamics wrt parameters.
%
% g = gpReversibleDynamicsLogLikeGradients(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpReversibleDynamicsLogLikeGradients.m version 1.1



if model.k ==0 & ~model.learn & ~model.learnScales
  g = [];
  return
end

g = gpLogLikeGradients(model);

if ~model.learn
  % If we aren't learning model parameters extract only X_u;
  % this is inefficient (but neater in the code) as we have also computed parameters 
  if ~model.learnScales
    g = g(1:model.k*model.q);
  else
    switch model.approx
     case 'ftc'
      g =  [g(end-model.d + 1:end)];
     case {'dtc', 'fitc', 'pitc'}
      g =  [g(1:model.k*model.q) g(end-model.d:end-1)];
    end
  end
end
