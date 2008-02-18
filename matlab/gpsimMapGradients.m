function g = gpsimMapGradients(param, model)

% GPSIMMAPGRADIENTS Compute the gradients of the log likelihood of a GPSIMMAP model.
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood,
% gpsimMapGradient, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% GPSIM  
  
  
g = zeros(1,length(param));
Nrep = length(model.comp);
   
for rep = 1:Nrep   % Work out likelihood gradient for each replicate
  options = defaultOptions;
  model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, param);
  model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);

  dg{rep} = gpsimMapLogLikeGradients(model.comp{rep});
  g = g - dg{rep};
end


