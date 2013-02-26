function kern = whiteblockKernParamInit(kern)

% WHITEBLOCKKERNPARAMINIT WHITE BLOCK kernel parameter initialisation.
%
%	Description:
%	The white block noise kernel arises from assuming independent Gaussian
%	noise for each point in the function for multiple outputs in which the
%	kernel matrix is computed in one unique block.
%	This kernel is not intended to be used independently, it is provided
%	so that it may be combined with the LMC kernel in a compound kernel.
%
%	KERN = WHITEBLOCKKERNPARAMINIT(KERN) initialises the white noise
%	block kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2010 Mauricio A. Alvarez


if isfield(kern, 'options') && isfield(kern.options, 'nout')
    kern.nout = kern.options.nout;
else
    % If the number of outputs is not provided, we assume is one
    kern.nout = 1;
end

kern.variance = exp(-2)*ones(1, kern.nout);
kern.nParams = kern.nout;

kern.transforms.index = 1:kern.nout;
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = true;
