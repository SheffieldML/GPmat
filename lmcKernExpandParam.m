function kern = lmcKernExpandParam(kern, params)

% LMCKERNEXPANDPARAM Expands parameters into a LMC kernel structure.
% FORMAT
% DESC returns a linear model of corregionalization kernel structure filled 
% with the parameters in the given vector. This is used as a helper 
% function to enable parameters to be optimised in, for example, the NETLAB
% optimisation functions. 
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

fhandle = str2func([kern.basicKernelType 'KernExpandParam']);
kern = fhandle(kern, params(1:kern.nParamsBK));
kern.A = reshape(params(kern.nParamsBK+1:end), kern.nout, kern.rankCorregMatrix);
kern.B = kern.A*kern.A';
