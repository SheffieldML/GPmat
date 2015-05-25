function marginalLL = gpsimMapMarginalLikelihood(model)

% GPSIMMAPMARGINALLIKELIHOOD The marginal likelihood approximation.
%
% COPYRIGHT : Pei Gao, 2008
  
% GPSIM

%/~
% ll = gpsimMapLogLikelihood(model);
% invC = model.invK*model.covf;
% marginalLL = ll - 0.5*model.f'*model.invK*model.f - 0.5*(model.logDetK ...
%      - model.logDetCovf);
%~/
marginalLL = log(det(eye(size(model.K))+model.K*model.W));
%/~
%  marginalLL = (model.logDetK - model.logDetCovf);
%   I = eye(size(model.K));
%   marginalLL = log(det(I+model.K*model.W));
%~/