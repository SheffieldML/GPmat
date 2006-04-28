function g = gpScaleBiasGradient(model)

% GPSCALEBIASGRADIENT Compute the gradient of the scale and bias.
%
% g = gpScaleBiasGradient(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpScaleBiasGradient.m version 1.1


g = [];
if model.learnScales
  g = 1./model.scale.*(model.innerProducts-1);
  fhandle = str2func([model.scaleTransform 'Transform']);
  g = g.*fhandle(model.scale, 'gradfact');
end
