function g = lfmLogLikeGradients(model)

% LFMLOGLIKEGRADIENTS Compute the gradients of the log likelihood of a LFM model.
% FORMAT
% DESC computes the gradients of the log likelihood of the given
% Gaussian process for use in a single input motif protein network.
% ARG model : the model for which the log likelihood is computed.
% RETURN g : the gradients of the parameters of the model.
% 
% SEEALSO : lfmCreate, lfmLogLikelihood, lfmGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
covGrad = 0.5*covGrad;
g = kernGradient(model.kern, model.t, covGrad);

%/~ In case we need priors in.
% Add contribution of any priors 
% g = g + kernPriorGradient(model.kern);
%~/


if isfield(model, 'fix')
  for i = 1:length(model.fix)
    g(model.fix(i).index) = 0;
  end
end
