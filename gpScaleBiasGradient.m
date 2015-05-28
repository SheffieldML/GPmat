function g = gpScaleBiasGradient(model)

% GPSCALEBIASGRADIENT Compute the log likelihood gradient wrt the scales.
% FORMAT
% DESC computes the gradient of the log likelihood with respect to the
% scales. In the future the gradients with respect to the biases
% may also be included. 
% ARG model : the model for which the gradients are to be computed.
% ARG g : the gradients of the likelihood with respect to the
% scales.
%
% SEEALSO : gpCreate, gpLogLikeGradients, gpLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP

g = [];
if model.learnScales
  g = 1./model.scale.*(model.innerProducts-1);
  fhandle = str2func([model.scaleTransform 'Transform']);
  g = g.*fhandle(model.scale, 'gradfact');
end
