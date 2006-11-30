function ll = gpsimMapFunctionalLogLikelihood(model)

% GPSIMMAPFUNCTIONALLOGLIKELIHOOD Compute the log likelihood of a GPSIMMAP model.
% FORMAT
% DESC computes the log likelihood of the given Gaussian process
% for use in a single input motif protein network. Acts as a
% wrapper for gpsimMapLogLikelihood.
% ARG model : the model for which the log likelihood is computed.
% RETURN ll : the log likelihood of the data set.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM

ll = gpsimMapLogLikelihood(model);