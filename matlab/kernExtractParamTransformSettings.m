function [paramtransformsettings, names] = kernExtractParamTransformSettings(kern)

% KERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings from kernel structure.
%
% FORMAT
%
% DESC Extract parameters from the kernel into a vector of
% parameters for optimisation. If the field 'transforms' is not 
% empty in the kernel, the parameters will be transformed back 
% into an unrestricted space suitable for optimisation. For
% example, if some parameters are required to be positive in the 
% kernel, and they are made positive by an exponential 
% transformation, then the parameters will be inverse-transformed 
% by applying the log function before returning. If any custom 
% settings have been provided for the transformations (like 
% custom ranges allowed for kernel parameters), they will be 
% taken into account in the inverse transformations.
%
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% 
% RETURN param : vector of parameters extracted from the kernel. 
%
% SEEALSO : kernExpandParam, scg, conjgrad
% 
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011.
%
% COPYRIGHT : Antti Honkela, 2012

% KERN

fname = [kern.type 'KernExtractParamTransformSettings'];

% Check if a specific method exists
if exist(fname),
  fhandle = str2func(fname);
  names = cell(1, kern.nParams);

  if nargout > 1
    [paramtransformsettings, names] = fhandle(kern);
  else
    paramtransformsettings = fhandle(kern);
  end
else % No, so let's improvise
  if nargout > 1,
    phandle = str2func([kern.type 'KernExtractParam']);
    [~, names] = phandle(kern);
  end

  [paramtransformsettings{1:length(kern.transforms)}] = ...
      deal(kern.transforms.transformsettings);
end
