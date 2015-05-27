function ll = gpsimMapLogLikelihood(model, df)

% GPSIMMAPLOGLIKELIHOOD Compute the log likelihood of a GPSIMMAP model.
% FORMAT
% DESC computes the log likelihood of the given Gaussian process
% for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN ll : the log likelihood of the data set.
% 
% SEEALSO : gpsimMapCreate, gpsimMapLogLikeGradient, gpsimMapObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML

numData = length(model.t);

% ll = model.f'*model.invK*model.f ...
%       + log(det(eye(size(model.W))+model.K*model.W));
ll = model.f'*model.invK*model.f ...
    - model.logDetCovf + model.logDetK;

if isfield(model, 'includeNoise') && model.includeNoise
  noiseMat = ones(length(model.t), 1)*model.noiseVar;
  yvar = model.yvar + noiseMat;
else
  yvar = model.yvar;
end

for i = 1:numData
  for j = 1:model.numGenes
    ind = i + (j-1)*numData;
    beta_ij = 1/yvar(ind);
    factor = (model.ypred(model.times_index(i), j)...
              - model.y(ind));
    ll = ll + factor*factor*beta_ij - log(beta_ij);
  end
end

ll = ll + numData*model.numGenes*log(2*pi);

ll = -.5*ll;
