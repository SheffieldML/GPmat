function model = gpReversibleDynamicsExpandParam(model, params)

% GPREVERSIBLEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
%
% model = gpReversibleDynamicsExpandParam(model, params)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpReversibleDynamicsExpandParam.m version 1.1



model = gpDynamicsExpandParam(model, params);
