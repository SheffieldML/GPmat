function model = robThreeDynamicsExpandParam(model, params)

% ROBTHREEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.

% SHEFFIELDML

if length(params)>0
  error(['There should be no placing of parameters in Robot Three ' ...
         'Dynamics']);
end
return;