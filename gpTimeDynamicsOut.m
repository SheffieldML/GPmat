function y = gpTimeDynamicsOut(model, x);

% GPTIMEDYNAMICSOUT Evaluate the output of GPTIMEDYNAMICS.
% FORMAT
% DESC evaluates the output of a given Gaussian process regressive dynamics model.
% ARG model : the model for which the output is being evaluated.
% ARG t : the time output is required.
% RETURN y : the output of the GP model. The function checks if
% there is a noise model associated with the GP, if there is, it is
% used, otherwise the mean of gpPosteriorMeanVar is returned.
%
% SEEALSO : gpTimeDynamicsCreate, gpOut
%
% COPYRIGHT : Neil D. Lawrence 2007

% FGPLVM

y = gpOut(model, x);
