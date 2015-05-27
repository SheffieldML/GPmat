function kern = indexardKernExpandParam(kern, params)

% INDEXARDKERNEXPANDPARAM Create kernel structure from INDEXARD kernel's parameters.
% FORMAT
% DESC returns a index based covariance function kernel structure filled with the
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
% SEEALSO : indexardKernParamInit, indexardKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

  kern.indexScales = params(1:end);
end
