function model = robTwoDynamicsExpandParam(model, params)

% ROBTWODYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
%
% model = robTwoDynamicsExpandParam(model, params)
%

% Copyright (c) 2006 Neil D. Lawrence
% robTwoDynamicsExpandParam.m version 1.1



if length(params)>0
  error(['There should be no placing of parameters in Robot Two ' ...
         'Dynamics']);
end
return;