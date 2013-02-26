function kern = gibbsKernSetLengthScaleFunc(kern, model)

% GIBBSKERNSETLENGTHSCALEFUNC Set the length scale function of the GIBBS kernel.
%
%	Description:
%
%	KERN = GIBBSKERNSETLENGTHSCALEFUNC(KERN, MODEL)
%	 Returns:
%	  KERN - the kernel with the given length scale function.
%	 Arguments:
%	  KERN - the GIBBS kernel for which you want to change the length
%	   scale function.
%	  MODEL - the length scale function you wish to use.
%	
%
%	See also
%	MODELCREATE, GIBBSKERNPARAMINIT


%	Copyright (c) 2006 Neil D. Lawrence



if model.inputDim ~= kern.inputDimension
  error('Length scale function must have same input dimension as kernel.');
end

kern.lengthScaleFunc = model;
kern.nParams = 1+kern.lengthScaleFunc.numParams;

