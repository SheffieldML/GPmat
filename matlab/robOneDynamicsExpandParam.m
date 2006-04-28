function model = robOneDynamicsExpandParam(model, params)

% ROBONEDYNAMICSEXPANDPARAM Place the parameters vector into the model for first robot dynamics.
%
% model = robOneDynamicsExpandParam(model, params)
%

% Copyright (c) 2006 Neil D. Lawrence
% robOneDynamicsExpandParam.m version 1.1



if length(params)>0
  error(['There should be no placing of parameters in Robot One ' ...
         'Dynamics']);
end
return;