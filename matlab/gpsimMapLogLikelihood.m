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

% GPSIM

ll = model.f'*model.invK*model.f ...
     - model.logDetCovf + model.logDetK;
for i = 1:length(model.t)
  for j = 1:model.numGenes
    beta_ij = 1/model.yvar(i+(j-1)*length(model.t));
    factor = (model.ypred(model.times_index(i), j)...
              - model.y(i+(j-1)*length(model.t)));
    ll = ll + factor*factor*beta_ij + log(beta_ij);
  end
end

ll = ll + length(model.t)*model.numGenes*log(2*pi);

ll = -.5*ll;
