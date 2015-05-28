function [f, g] = gpsimMapGradFuncWrapper(param, model)
  
% GPSIMMAPGRADFUNCWRAPPER wraps the log-likelihood function and the gradient function 
% together for minimize optimisation.

% FORMAT
% ARG param : a vector of parameters from the model.
% ARG model : the model.
% RETURN f : the log likelihood of the data set.
% RETURN g : gradient function of the model with respect to the parameters.  
  
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood, gpsimMapLogLikeGradients
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008
  
% SHEFFIELDML

g = zeros(1,length(param));
if isfield(model, 'comp')
    Nrep = length(model.comp);
    for rep=1:Nrep  %Work out likelihood gradient for each replicate
        options = defaultOptions;
        model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, param);
        model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);
        ll(rep) = gpsimMapLogLikelihood(model.comp{rep});
        if nargout > 1
          dg{rep} = gpsimMapLogLikeGradients(model.comp{rep});
          g = g - dg{rep};
        end
        fprintf('log-likelihood %2.4f;\t', ll(rep));
    end
else
    Nrep = 1;
    options = defaultOptions;
    model = gpsimMapExpandParam(model, param);
    model = gpsimMapUpdateF(model, options);
    ll = gpsimMapLogLikelihood(model);
    dg = gpsimMapLogLikeGradients(model);
    g = -dg;
    fprintf('log-likelihood %2.4f;\t', ll);
end

fprintf('\n');
f = - sum(ll);

if nargout>1 
  g = g';
end


