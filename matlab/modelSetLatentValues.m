function model = modelSetLatentValues(model, varargin)

% MODELSETLATENTVALUES Set the latent variables for dynamics models in the GPLVM.
%
% model = modelSetLatentValues(model, varargin)
%

% Copyright (c) 2006 Neil D. Lawrence
% modelSetLatentValues.m version 1.1


fhandle = str2func([model.type 'SetLatentValues']);
model = fhandle(model, varargin{:});