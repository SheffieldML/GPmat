function model = modelSetLatentValues(model, varargin)

% MODELSETLATENTVALUES Set the latent variables for dynamics models in the GPLVM.
% FGPLVM

fhandle = str2func([model.type 'SetLatentValues']);
model = fhandle(model, varargin{:});
