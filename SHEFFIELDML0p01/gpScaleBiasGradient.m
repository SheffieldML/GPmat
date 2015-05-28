function g = gpScaleBiasGradient(model)

% GPSCALEBIASGRADIENT Compute the log likelihood gradient wrt the scales.
%
%	Description:
%
%	GPSCALEBIASGRADIENT(MODEL, G) computes the gradient of the log
%	likelihood with respect to the scales. In the future the gradients
%	with respect to the biases may also be included.
%	 Arguments:
%	  MODEL - the model for which the gradients are to be computed.
%	  G - the gradients of the likelihood with respect to the scales.
%	
%
%	See also
%	GPCREATE, GPLOGLIKEGRADIENTS, GPLOGLIKELIHOOD


%	Copyright (c) 2006 Neil D. Lawrence


g = [];
if model.learnScales
  g = 1./model.scale.*(model.innerProducts-1);
  fhandle = str2func([model.scaleTransform 'Transform']);
  g = g.*fhandle(model.scale, 'gradfact');
end
