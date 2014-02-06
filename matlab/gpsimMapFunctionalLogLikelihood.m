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
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

ll = gpsimMapLogLikelihood(model)-0.5*model.logDetCovf;

% Add constraints
if isfield(model, 'priorProtein') && ~isempty(model.priorProtein)
  nCons = length(model.priorProtein);
  for i=1:nCons
    ftimeIndex = find((model.priorProteinTimes(i)-model.mapt)==0);
    ll = ll - 0.5*model.consLambda*(model.f(ftimeIndex) - model.priorProtein(i))^2;
  end
end
