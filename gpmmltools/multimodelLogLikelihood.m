function ll = multimodelLogLikelihood(model)

% MULTIMODELLOGLIKELIHOOD Log likelihood of MULTIMODEL model.
% FORMAT
% DESC computes the log likelihood of  the multi-task learning wrapper model.
% ARG model : the model structure for which log likelihood is being computed.
% RETURN ll : the computed log likelihood.
%
% SEEALSO : multimodelCreate, multimodelLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

ll = 0;
for i = 1:length(model.comp)
  ll = ll + modelLogLikelihood(model.comp{i});
end
