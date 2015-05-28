function [samples, E_hist, g_mean] = gpnddisimSampleHMC(model, display, iters, options);

% GPNDDISIMSAMPLEHMC Do HMC sampling for the GPNDDISIM model.
% FORMAT
% DESC performs HMC sampling for the Gaussian process single input
% motif model for a given number of iterations.
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation
% proceeds, set to 2 for the most verbose and 0 for the least
% verbose.
% ARG iters : number of samples to return.
% RETURN samples : the samples.
%
% SEEALSO : hmc, gpsimCreate, gpsimGradient, gpsimObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2013
%
% COPYRIGHT : Jaakko Peltonen, 2011

% GPNDDISIM


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end


[params,paramnames] = modelExtractParam(model);

if nargin < 4,
  options = hmcDefaultOptions;
end

if isfield(model, 'uniformPriors') && model.uniformPriors,
  options.iters = iters;
  bounds = gpnddisimExtractParamTransformSettings(model);
  bounds = cat(1, bounds{:})';
  [samples, E_hist, g_mean] = ...
      myhmc_bounded('modelObjective', params,  options, ...
		    'modelGradient', bounds, model);
else
  hmcopt = optOptions;
  if display
    hmcopt(1) = display;
    if length(params) <= 100
      hmcopt(9) = 1;
    end
  end

  % Momentum persistence
  hmcopt(5) = 1;

  % Leapfrog steps
  hmcopt(7) = options.tau;

  % hmcopt(18) = epsilon, step length
  hmcopt(18) = options.epsilon;

  % Number of samples to return
  hmcopt(14) = iters;

  samples = hmc('modelObjective', params,  hmcopt, ...
		'modelGradient', model);

  E_hist = [];
  g_mean = [];
end
