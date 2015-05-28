function model = robTwoDynamicsExpandParam(model, params)

% ROBTWODYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
%
%	Description:
%	model = robTwoDynamicsExpandParam(model, params)
%

if length(params)>0
  error(['There should be no placing of parameters in Robot Two ' ...
         'Dynamics']);
end
return;