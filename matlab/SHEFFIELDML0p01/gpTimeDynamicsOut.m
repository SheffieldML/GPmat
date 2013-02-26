function y = gpTimeDynamicsOut(model, x);

% GPTIMEDYNAMICSOUT Evaluate the output of GPTIMEDYNAMICS.
%
%	Description:
%
%	Y = GPTIMEDYNAMICSOUT(MODEL, T) evaluates the output of a given
%	Gaussian process regressive dynamics model.
%	 Returns:
%	  Y - the output of the GP model. The function checks if there is a
%	   noise model associated with the GP, if there is, it is used,
%	   otherwise the mean of gpPosteriorMeanVar is returned.
%	 Arguments:
%	  MODEL - the model for which the output is being evaluated.
%	  T - the time output is required.
%	
%
%	See also
%	GPTIMEDYNAMICSCREATE, GPOUT


%	Copyright (c)  Neil D. Lawrence 2007


y = gpOut(model, x);