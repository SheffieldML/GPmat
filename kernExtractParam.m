function [params, names] = kernExtractParam(kern)

% KERNEXTRACTPARAM Extract parameters from kernel structure.
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

% KERN

fhandle = str2func([kern.type 'KernExtractParam']);
names = cell(1, kern.nParams);

if nargout > 1
  [params, names] = fhandle(kern);
else
  params = fhandle(kern);
end
  
%/~
if any(isnan(params));
  warning('Parameter has gone to NaN')
end
%~/

% Check if parameters are being optimised in a transformed space,
% and if they are, apply the inverse transformations to get the
% parameters back into an unrestricted space.

% untransformed_params=params;
if ~isempty(kern.transforms)
  % Process each transformation. Each transformation may affect
  % several parameters.
  for i = 1:length(kern.transforms)
    % Get the parameter indices affected by the i:th transformation.
    index = kern.transforms(i).index;
    
    % Get the function handle for the i:th transformation
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    
    % If custom settings have been provided for the i:th
    % transformation, use them, otherwise call the transformation
    % without parameters.
    if isfield( kern.transforms(i),'transformsettings' ) && ~isempty(kern.transforms(i).transformsettings')
      params(index) = fhandle(params(index), 'xtoa', kern.transforms(i).transformsettings);
    else
      params(index) = fhandle(params(index), 'xtoa');
    end    
  end
end

%fprintf(1, 'after transformation:\n');
%kern.type
%kern.transforms
%for k=1:length(kern.transforms),
%  kern.transforms(k).type
%  kern.transforms(k).index
%  kern.transforms(k).transformsettings
%  
%  l=kern.transforms(k).index;
%  fprintf(1, '%d: %f --> %f (%f %f)\n', l, untransformed_params(l), params(l), kern.transforms(k).transformsettings(1), kern.transforms(k).transformsettings(2));
%end;

