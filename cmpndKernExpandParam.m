function kern = cmpndKernExpandParam(kern, params)


% CMPNDKERNEXPANDPARAM Create kernel structure from CMPND kernel's parameters.
% FORMAT
% DESC returns a compound kernel structure filled with the
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
% SEEALSO : cmpndKernParamInit, cmpndKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN

params = params*kern.paramGroups';
startVal = 1;
endVal = 0;
kern.whiteVariance = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  kern.comp{i} = kernExpandParam(kern.comp{i}, params(1, startVal:endVal));
  startVal = endVal + 1;
  if strcmp(kern.comp{i}.type(1:min([end 5])), 'white')
    % If kernel name starts with white, assume it is a white noise term.
    kern.whiteVariance = kern.whiteVariance + kern.comp{i}.variance;
  else
    if(isfield(kern.comp{i}, 'whiteVariance'))
      kern.whiteVariance = kern.whiteVariance + ...
          kern.comp{i}.whiteVariance;
    end
  end
end
