function g = negNoiseGradientParam(params, model, prior)

% NEGNOISEGRADIENTPARAM Wrapper function for calling noise gradients.
%
%	Description:
%	g = negNoiseGradientParam(params, model, prior)
%

model.noise = noiseExpandParam(model.noise, params);
fhandle = str2func([model.noise.type 'GradientParam']);
g = - fhandle(model);

if prior
  g = g + params;
end