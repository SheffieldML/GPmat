function kern = lfmKernExpandParam(kern, params)

% LFMKERNEXPANDPARAM Create kernel structure from LFM kernel's parameters.
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
% SEEALSO : lfmKernParamInit, lfmKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% MODIFICATIONS : David Luengo, 2008, 2009

% KERN

kern.mass = params(1);
kern.spring = params(2);
kern.damper = params(3);
kern.inverseWidth = params(4);
kern.sensitivity = params(5);

kern.alpha = kern.damper./(2*kern.mass);
kern.omega = sqrt(kern.spring./kern.mass-kern.alpha.^2);
kern.gamma = kern.alpha + j*kern.omega;

kern.zeta = kern.damper./(2*sqrt(kern.mass*kern.spring));
kern.omega_0 = sqrt(kern.spring./kern.mass);
