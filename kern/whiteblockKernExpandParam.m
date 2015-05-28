function kern = whiteblockKernExpandParam(kern, params)

% WHITEBLOCKKERNEXPANDPARAM Fill WHITEBLOCK kernel structure with params.
% FORMAT
% DESC returns a white noise block kernel structure filled with the
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
% SEEALSO : whiteblockKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

kern.variance = params;
