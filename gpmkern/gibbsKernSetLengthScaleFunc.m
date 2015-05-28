function kern = gibbsKernSetLengthScaleFunc(kern, model)

% GIBBSKERNSETLENGTHSCALEFUNC Set the length scale function of the GIBBS kernel.
% FORMAT
% ARG kern : the GIBBS kernel for which you want to change the
% length scale function.
% ARG model : the length scale function you wish to use.
% RETURN kern : the kernel with the given length scale function.
%
% SEEALSO : modelCreate, gibbsKernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN


if model.inputDim ~= kern.inputDimension
  error('Length scale function must have same input dimension as kernel.');
end

kern.lengthScaleFunc = model;
kern.nParams = 1+kern.lengthScaleFunc.numParams;

