function ll = pmvuLogLikelihood(model)

% PMVULOGLIKELIHOOD Log likelihood of PMVU model.
% FORMAT
% DESC computes the log likelihood of  the probabilistic maximum variance unfolding model.
% ARG model : the model structure for which log likelihood is being computed.
% RETURN ll : the computed log likelihood.
%
% SEEALSO : pmvuCreate, pmvuLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence 2009
 
% MLTOOLS

ll = -model.d*model.logDetK;
ll = ll - sum(sum(model.Y.*(model.L*model.Y))) -model.gamma*model.traceY;
ll = ll*0.5;