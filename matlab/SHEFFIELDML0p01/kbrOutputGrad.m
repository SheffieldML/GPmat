function g = kbrOutputGrad(model, X, dim)

% KBROUTPUTGRAD Evaluate derivatives of KBR model outputs with respect to parameters.
%
%	Description:
%
%	G = KBROUTPUTGRAD(MODEL, X) evaluates the derivates of a kernel
%	based regression outputs with respect to the parameters of the
%	multi-layer perceptron.
%	 Returns:
%	  G - the gradient of the outputs of the kernel based regression
%	   with respect to each of the parameters. The size of the matrix is
%	   number of data x number of parameters x number of outputs of the
%	   model.
%	 Arguments:
%	  MODEL - the model for which the derivatives are to be computed.
%	  X - the input data locations where the gradients are to be
%	   computed.
%	
%	
%
%	See also
%	KBRCREATE, KBRDERIV


%	Copyright (c) 2008 Neil D. Lawrence


%	With modifications by Carl Henrik Ek 2008


numData = size(X, 1);
if(nargin<=2)
  for i = 1:model.outputDim
    startZeros = zeros(numData, numData*(i - 1));
    finishZeros = zeros(numData, numData*(model.outputDim-i));
    startZeros2 = zeros(numData, (i - 1));
    finishZeros2 = zeros(numData, (model.outputDim-i));
    g(:, :, i) = [startZeros model.K finishZeros startZeros2 ones(numData, 1) finishZeros2];
  end
else
  g(:,:) = [model.K ones(numData,1)];
end