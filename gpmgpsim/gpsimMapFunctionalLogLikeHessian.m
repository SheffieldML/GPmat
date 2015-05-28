function H = gpsimMapFunctionalLogLikeHessian(model)

% GPSIMMAPFUNCTIONALLOGLIKEHESSIAN Compute the functional Hessian for GPSIMMAP.
% FORMAT
% DESC computes the functional Hessian of the log likelihood for
% use in the MAP approximation to the GPSIM posterior solution.
% ARG model : the model for which the Hessian is to be computed.
% RETURN H : the Hessian of the log likelihood with respect to the
% points of the function.
% 
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalUpdateW, gpsimMapUpdateKernels
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006

% SHEFFIELDML

H = - model.W - model.invK;
