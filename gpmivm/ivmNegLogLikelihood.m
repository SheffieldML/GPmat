function e = ivmNegLogLikelihood(params, model)

% IVMNEGLOGLIKELIHOOD Wrapper function for calling IVM likelihood.
% FORMAT
% DESC is a wrapper function for calling the IVM log likelihood.
% ARG param : the parameters where the log likelihood is to be
% evaluated.
% ARG model : the model structure for which the log likelihood is
% being evaluated.
% ARG e : the negative log likelihood of the model.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : noiseExpandParam, ivmLogLikelihood

% IVM

model.noise = noiseExpandParam(model.noise, params);
e = - ivmLogLikelihood(model);
