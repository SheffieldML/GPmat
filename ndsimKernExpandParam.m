function kern = ndsimKernExpandParam(kern, params)

% NDSIMKERNEXPANDPARAM Create kernel structure from NDSIM kernel's parameters.
% FORMAT
% DESC returns a single input motif kernel structure with zero
% decay, filled with the
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
% SEEALSO : ndsimKernParamInit, ndsimKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

kern.inverseWidth = params(1);
%fprintf(1,'NDSIM kern.inverseWidth=%f\n',kern.inverseWidth);

if isfield(kern, 'isNegativeS') && kern.isNegativeS
  kern.sensitivity = params(2);
  kern.variance = kern.sensitivity*kern.sensitivity;
else
  kern.variance = params(2);
end
if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  kern.initialVariance = params(3);
end
