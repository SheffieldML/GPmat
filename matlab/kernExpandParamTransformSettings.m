function kern = kernExpandParamTransformSettings(kern, paramtransformsettings)

% KERNEXPANDPARAMTRANSFORMSETTINGS Expand parameters' transform settings to form a kernel structure.
%
% FORMAT
%
% DESC returns a kernel structure filled with the parameter
% transform settings in the given cell array. This is used as a 
% helper function to enable parameters to be optimised in, for 
% example, the NETLAB optimisation functions.
%
% ARG kern : the kernel structure in which the parameter
% transformation settings are to be placed.
%
% ARG paramtransformsettings : cell array (nParams x 1) of
% parameter transform settings which are to be placed in the
% kernel structure. Each transformation setting can be of any
% form, depending on what information the corresponding
% transformation needs. For example, a transformation setting
% might be a desired output range like [0 100], or it could
% be some more complicated structure.
%
% RETURN kern : kernel structure with the given parameter 
% transform settings in the relevant locations.
%
% SEEALSO : kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% COPYRIGHT : Antti Honkela, 2011

% KERN

%fprintf(1,'kernExpandParamTransformSettings step1\n');

% If no transformation settings are provided, do nothing, otherwise
% simply call the kernel-specific function.
if ~isempty(paramtransformsettings),
  % Get the handle of the kernel-specific function and call it
  fhandle = str2func([kern.type 'KernExpandParamTransformSettings']);
  kern = fhandle(kern, paramtransformsettings);
end;

if isfield(kern, 'priors')
  if iscell(kern.priors)
    for k=1:length(kern.priors),
      if isfield(kern.priors{k}, 'isBounded') && ...
	    kern.priors{k}.isBounded,
	l = kern.priors{k}.index;
	if length(l) > 1,
	  error('Unable to set bounds to a prior shared by multiple params')
	end
	kern.priors{k} = priorSetBounds(kern.priors{k}, paramtransformsettings{l});
      end
    end
  else
    for k=1:length(kern.priors),
      if isfield(kern.priors(k), 'isBounded') && ...
	    kern.priors(k).isBounded,
	l = kern.priors(k).index;
	if length(l) > 1,
	  error('Unable to set bounds to a prior shared by multiple params')
	end
	kern.priors(k) = priorSetBounds(kern.priors(k), paramtransformsettings{l});
      end
    end
  end
end

%fprintf(1,'kernExpandParamTransformSettings done\n');
