function g = gpnddisimGradient(params, model)

% GPDISIMGRADIENT Gradient wrapper for a GPDISIM model.
% FORMAT 
% DESC wraps the log likelihood gradient function to return the
% gradient of the negative of the log likelihood. This can then be
% used in, for example, NETLAB, minimisation tools.
% ARG params : the parameters of the model.
% ARG model : the model for which gradients will be computed.
% RETURN g : the returned gradient of the negative log likelihood
% for the given parameters.
%
% SEEALSO : scg, conjgrad, gpsimCreate, gpsimObjective, gpsimLogLikeGradient, gpsimOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007
  
% GPSIM

try,
    model = gpnddisimExpandParam(model, params);
    g = - gpnddisimLogLikeGradients(model);
catch,
    g = 0*gpnddisimLogLikeGradients(model);
end

%plotpredictions(model,[0:5:1280]',2,1,1,'exampletitle');
%drawnow;

%g


% g=-testModelGradientBasic(model,1,1);
