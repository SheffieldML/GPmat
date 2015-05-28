function samples = gpsimSampleHMC(model, display, iters);

% GPSIMSAMPLEHMC Do HMC sampling for the GPSIM model.
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
% MODIFICATIONS : Antti Honkela, 2007

% SHEFFIELDML


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end


params = modelExtractParam(model);

options = optOptions;
if display
  options(1) = display;
  if length(params) <= 100
    options(9) = 1;
  end
end
% Momentum persistence
options(5) = 1;
% Leapfrog steps
options(7) = 10;
options(18) = .01;
% Number of samples to return
options(14) = iters;

samples = hmc('modelObjective', params,  options, ...
	      'modelGradient', model);
