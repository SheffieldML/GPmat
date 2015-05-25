function kern = simKernExpandParam(kern, params)

% SIMKERNEXPANDPARAM Create kernel structure from SIM kernel's parameters.
% FORMAT
% DESC returns a single input motif kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : simKernParamInit, simKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.decay = params(1);
kern.inverseWidth = params(2);
if isfield(kern, 'isNegativeS') && kern.isNegativeS
  kern.sensitivity = params(3);
  kern.variance = kern.sensitivity*kern.sensitivity;
else
  kern.variance = params(3);
end
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  kern.initialVariance = params(4);
end
