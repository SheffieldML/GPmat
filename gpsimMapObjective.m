function f = gpsimMapObjective(param, model)
  
% GPSIMMAPOBJECTIVE Compute the objective function  of the GPSIMMAP model.
% FORMAT
% DESC computes the objective function of the model given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the objective function is computed.
% RETURN f : the objective function of the model.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikelihood,
% gpsimMapGradient, gpsimMapFunctionalLogLikeGradients
%
% COPYRIGHT : Pei Gao, Magnus Rattray and Neil D. Lawrence, 2008

% GPSIM   
  

if isfield(model, 'comp')
    Nrep = length(model.comp);
    for rep=1:Nrep  %Work out likelihood gradient for each replicate
        options = defaultOptions;
        model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, param);
%         f = gpsimMapFunctionalExtractParam(model.comp{rep});
%         model.comp{rep} = gpsimMapFunctionalExpandParam(model.comp{rep}, f);        
%        options(1) = 1;
        model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);
        ll(rep) = gpsimMapLogLikelihood(model.comp{rep});
    end
end

f = - sum(ll);

