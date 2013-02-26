function kern = multiKernExpandParam(kern, params)

% MULTIKERNEXPANDPARAM Create kernel structure from MULTI kernel's parameters.
%
%	Description:
%
%	KERN = MULTIKERNEXPANDPARAM(KERN, PARAM) returns a multiple output
%	block kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
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
%	MULTIKERNPARAMINIT, MULTIKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2006 Neil D. Lawrence


params = params*kern.paramGroups';
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isfield(kern, 'fixedBlocks') || ~kern.fixedBlocks(i),
    kern.comp{i} = kernExpandParam(kern.comp{i}, params(1, startVal:endVal));
  end
  startVal = endVal + 1;
end


