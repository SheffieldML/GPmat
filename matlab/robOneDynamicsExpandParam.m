function model = robOneDynamicsExpandParam(model, params)

% ROBONEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.

% FGPLVM

if length(params)>0
  error(['There should be no placing of parameters in Robot One ' ...
         'Dynamics']);
end
return;
