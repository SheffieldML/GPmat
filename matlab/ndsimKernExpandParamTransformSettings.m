function kern = ndsimKernExpandParamTransformSettings(kern, paramtransformsettings)

% NDSIMKERNEXPANDPARAMTRANSFORMSETTINGS Create kernel structure from SIM kernel's parameters' transform settings.
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
% SEEALSO : simKernParamInit, simKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
% 
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN


%fprintf(1,'expanding SIM kernel param-transforms settings:\n');
%paramtransformsettings
%for k=1:length(paramtransformsettings),
%  paramtransformsettings{k}
%end;
%pause

if isfield(kern, 'options') ...
      && isfield(kern.options, 'isNegativeS') ...
      && kern.options.isNegativeS,  
  if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
    % positive transformations for: inversewidth, initial variance
    for k=1:2,
      kern.transforms(k).transformsettings = paramtransformsettings{k};
    end;
  else
    % positive transformations for: inversewidth
    for k=1:1,
      kern.transforms(k).transformsettings = paramtransformsettings{k};
    end;
  end;
    
else
  if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
    % positive transformations for: inversewidth, kernel-variance,
    % initial variance
    for k=1:3,
      kern.transforms(k).transformsettings = paramtransformsettings{k};
    end;
  else
    % positive transformations for: inversewidth, kernel-variance
    for k=1:2,
      kern.transforms(k).transformsettings = paramtransformsettings{k};
    end;
  end;
end;

