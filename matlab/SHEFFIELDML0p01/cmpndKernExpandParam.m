function kern = cmpndKernExpandParam(kern, params)

% CMPNDKERNEXPANDPARAM Create kernel structure from CMPND kernel's parameters.
%
%	Description:
%
%	KERN = CMPNDKERNEXPANDPARAM(KERN, PARAM) returns a compound kernel
%	structure filled with the parameters in the given vector. This is
%	used as a helper function to enable parameters to be optimised in,
%	for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	   placed.
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%
%	See also
%	CMPNDKERNPARAMINIT, CMPNDKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


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
