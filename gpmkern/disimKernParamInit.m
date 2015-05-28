function kern = disimKernParamInit(kern)

% DISIMKERNPARAMINIT DISIM kernel parameter initialisation.
% The driven input single input motif (DISIM) kernel is specifically designed for
% working with gene networks where there is assumed to be a single
% transcription factor controlling several genes. This transcription
% factor, in turn, is driven by its own gene's RNA. The model takes the
% following form: each gene is
% related to the transcription factor through the following
% differential equation,
%
% dx(t)/dt = B + C f(t-delta) - D x(t),
%
% where D is a decay term, C is a response term, delta is a time delay
% and B is an initial level. Then if f(t) is assumed to be the result of
% a further differential equation,
% 
% df(t)/dt = Sx'(t) - D' x'(t) 
%
% where x'(t) is assumed to come from a Gaussian process with an RBF
% covariance function f(t) is a Gaussian process with a covariance function
% provided by the single input motif kernel (SIM) and x(t) is a Gaussian
% process with covariance function provided by this kernel, the DISIM kernel.
%
% The kernel is designed to interoperate with the multiple output
% block kernel so that f(t) can be inferred given several different
% instantiations of x(t) (associated with different genes).
%
% The parameters (B, C, delta, S, D and D') are constrained positive.
%
% FORMAT
% DESC initialises the single input motif
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit, simKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007, 2009

% KERN

if kern.inputDimension > 1
  error('DISIM kernel only valid for one-D input.')
end

if isfield(kern, 'options') && isfield(kern.options, 'gaussianInitial') && ...
      kern.options.gaussianInitial,
  kern.gaussianInitial = 1;
  kern.initialVariance = 1;
else
  kern.gaussianInitial = 0;
end

kern.di_decay = .1;
kern.inverseWidth = 1;
kern.di_variance = 1;
kern.decay = 1;
kern.variance = 1;
kern.rbf_variance = 1;

if kern.gaussianInitial,
  kern.nParams = 7;
else
  kern.nParams = 6;
end

kern.transforms.index = 1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
