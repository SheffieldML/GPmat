function kern = heatKernExpandParam(kern, params)

% HEATKERNEXPANDPARAM Create kernel structure from HEAT kernel's parameters.
% FORMAT
% DESC returns a heat kernel structure filled with the parameters in the 
% given vector. This is used as a helper function to enable parameters to 
% be optimised in, for example, the NETLAB optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : heatKernParamInit, heatKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

kern.decay = params(1);
kern.diffusion = params(2);
kern.inverseWidthTime = params(3);
kern.inverseWidthSpace = params(4);
kern.sensitivity = params(5);
kern.sim.inverseWidth = kern.inverseWidthTime;
if kern.includeIC
   kern.inverseWidthSpaceIC = params(6);
   kern.sensitivityIC = params(7);
end

% if kern.includeIndSens
%     if kern.includeIC
%         kern.sensitivitySpace = params(8);
%     else
%         kern.sensitivitySpace = params(6);
%     end
% end
