function model = robThreeDynamicsExpandParam(model, params)

% ROBTHREEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
%
% model = robThreeDynamicsExpandParam(model, params)
%

% Copyright (c) 2006 Neil D. Lawrence
% robThreeDynamicsExpandParam.m version 



if length(params)>0
  error(['There should be no placing of parameters in Robot Three ' ...
         'Dynamics']);
end
return;