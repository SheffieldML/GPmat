function factors = kernFactors(kern, factorType)

% KERNFACTORS Extract factors associated with transformed
% optimisation space. "factorType" is one of 'atox', 'xtoa',
% or 'gradfact', as follows.
%
% 'atox': transform unrestricted input parameters to a suitably
% restricted output range.
%
% 'xtoa': transform restricted parameters back to the unrestricted
% space where e.g. gradient-based parameter optimization is done.
%
% 'gradfact': These factors are derivatives of the 
% form dx/da, where "x" is the transformed parameter (for example
% a parameter that is restricted to be positive by a suitable
% transformation) and "a" is the untransformed parameter, which
% usually can freely take any real values.
%
% COPYRIGHT : unknown original copyright, possibly Neil Lawrence
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

factors.index = [];
factors.val = [];
if ~isempty(kern.transforms)
  fhandle = str2func([kern.type 'KernExtractParam']);
  params = fhandle(kern);
  
  % Process each transformation used with this kernel. Each
  % transformation may affect several parameters.
  for i = 1:length(kern.transforms)
    % Get the parameter indices involved in the i:th transformation
    index = kern.transforms(i).index;
    factors.index = [factors.index index];
    
    % Get the handle of the transformation function of the i:th transformation
    fhandle = str2func([kern.transforms(i).type 'Transform']);
    
    % If the transformation has been provided with specific
    % settings (such as a custom output range), use the settings, 
    % otherwise transform without settings
    if isfield(kern.transforms(i),'transformsettings') && ~isempty(kern.transforms(i).transformsettings)
      factors.val = [factors.val  ...
		     fhandle(params(index), factorType, kern.transforms(i).transformsettings)];
    else
      factors.val = [factors.val  ...
		     fhandle(params(index), factorType)];
    end
    
  end
end
