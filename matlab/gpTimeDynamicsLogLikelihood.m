function ll = gpTimeDynamicsLogLikelihood(model)

% GPTIMEDYNAMICSLOGLIKELIHOOD Give the log likelihood of GP time dynamics.
% FORMAT
% DESC Computes the log likelihood of GP time dynamics in a GP-LVM model.
% ARG model : the GP model for which log likelihood is to be
% computed.
% RETURN ll : the log likelihood of the data in the GP model.
%
% SEEALSO : gpTimeDynamicsCreate, gpTimeDynamicsLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

ll = gpLogLikelihood(model);

%/~
% % Use prior on first point in each sequence only for dynamics models.
% if isfield(model, 'prior') &  ~isempty(model.prior)
%   ll = ll + priorLogProb(model.prior, model.X(1, :));
%   for i = 2:length(model.seq)
%     ll = ll + priorLogProb(model.prior, model.X(model.seq(i-1)-i,:));
%   end
% end  
%~/