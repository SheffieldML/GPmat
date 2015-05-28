function ll = linearLogLikelihood(model)

% LINEARLOGLIKELIHOOD Linear model log likelihood.
% FORMAT
% DESC computes the log likelihood of a linear model. 
% ARG model : the model structure for computing the log likelihood.
% RETURN ll : the model log likelihood.
%
% SEEALSO : modelLogLikeihood, mlperr
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

N = size(model.y, 1);
centred = model.y - repmat(model.b, N, 1) - model.X* ...
       model.W;

centred = centred*model.beta;
ll = N*model.outputDim*(log(model.beta) - log(2*pi)) - sum(sum(centred.*centred))* ...
     model.beta;
ll = ll*0.5;
