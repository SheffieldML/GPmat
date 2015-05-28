function kern = fileKernParamInit(kern)

% FILEKERNPARAMINIT FILE kernel parameter initialisation.
% The file (FILE) kernel is designed for working with pre-computed
% kernel files that are saved on disk. It loads in a file the
% first time it is accessed and then caches it in memory (in single
% precision). The cacheing is done in by the fileKernRead command.
%
% This kernel type was originally written for working with Bill
% Stafford Noble's yeast data, available at
% http://noble.gs.washington.edu/proj/yeast/
%
% SEEALSO : fileKernRead
%
% FORMAT
% DESC initialises the stored file
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = [1];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
