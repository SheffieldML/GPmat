function ll = gpDynamicsLogLikelihood(model)

% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of GP dynamics.
%
%	Description:
%
%	LL = GPDYNAMICSLOGLIKELIHOOD(MODEL) Computes the log likelihood of
%	GP dynamics in a GP-LVM model.
%	 Returns:
%	  LL - the log likelihood of the data in the GP model.
%	 Arguments:
%	  MODEL - the GP model for which log likelihood is to be computed.
%	
%	
%	
%
%	See also
%	GPDYNAMICSCREATE, GPDYNAMICSLOGLIKEGRADIENTS, MODELLOGLIKELIHOOD


%	Copyright (c) 2006 Neil D. Lawrence and Carl Henrik Ek
%	Copyright (c) 2008 Carl Henrik Ek


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
