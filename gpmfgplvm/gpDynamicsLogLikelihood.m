function ll = gpDynamicsLogLikelihood(model)

% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of GP dynamics.
% FORMAT
% DESC Computes the log likelihood of GP dynamics in a GP-LVM model.
% ARG model : the GP model for which log likelihood is to be
% computed.
% RETURN ll : the log likelihood of the data in the GP model.
%
% SEEALSO : gpDynamicsCreate, gpDynamicsLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence and Carl Henrik Ek, 2006
%
% COPYRIGHT : Carl Henrik Ek, 2008
%

% FGPLVM

ll = gpLogLikelihood(model);

% Use prior on first point in each sequence only for dynamics models.
if isfield(model, 'prior') &  ~isempty(model.prior)
  if(isfield(model,'indexIn')&&~isempty(model.indexIn))
    ll = ll + priorLogProb(model.prior, model.X(1,model.indexIn));
     for i = 2:length(model.seq)
       ll = ll + priorLogProb(model.prior, model.X(model.seq(i-1)-i,model.indexIn));
     end
  else
    ll = ll + priorLogProb(model.prior, model.X(1, :));
    for i = 2:length(model.seq)
      ll = ll + priorLogProb(model.prior, model.X(model.seq(i-1)-i,:));
    end
  end
end  
