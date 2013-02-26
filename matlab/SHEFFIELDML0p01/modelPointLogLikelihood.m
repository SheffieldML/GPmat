function ll = modelPointLogLikelihood(model, varargin)

% MODELPOINTLOGLIKELIHOOD Compute the log likelihood of a given point.
%
%	Description:
%
%	LL = MODELPOINTLOGLIKELIHOOD(MODEL, ...) computes the log likelihood
%	of the given model.
%	 Returns:
%	  LL - the log likelihood of the given data point.
%	 Arguments:
%	  MODEL - the model for which the log likelihood is to be computed.
%	  ... - additional arguments as required.
%	
%
%	See also
%	MODELLOGLIKELIHOOD, MODELCREATE


%	Copyright (c) 2006 Neil D. Lawrence


fhandle = str2func([model.type 'PointLogLikelihood']);
ll = fhandle(model, varargin{:});