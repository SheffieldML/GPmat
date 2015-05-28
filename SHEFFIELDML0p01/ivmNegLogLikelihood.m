function e = ivmNegLogLikelihood(params, model)

% IVMNEGLOGLIKELIHOOD Wrapper function for calling IVM likelihood.
%
%	Description:
%
%	IVMNEGLOGLIKELIHOOD(PARAM, MODEL, E) is a wrapper function for
%	calling the IVM log likelihood.
%	 Arguments:
%	  PARAM - the parameters where the log likelihood is to be
%	   evaluated.
%	  MODEL - the model structure for which the log likelihood is being
%	   evaluated.
%	  E - the negative log likelihood of the model.
%	
%
%	See also
%	NOISEEXPANDPARAM, IVMLOGLIKELIHOOD


%	Copyright (c) 2005 Neil D. Lawrence


model.noise = noiseExpandParam(model.noise, params);
e = - ivmLogLikelihood(model);
