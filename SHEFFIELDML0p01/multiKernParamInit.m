function kern = multiKernParamInit(kern)

% MULTIKERNPARAMINIT MULTI kernel parameter initialisation.
%
%	Description:
%	The multiple output block kernel (MULTI) is a wrapper kernel
%	designed to represent the situation where there are several Gaussian
%	processes with correlated outputs. The kernel stores the parameters
%	of the kernels which form the block diagonal of the covariance
%	matrix of the full process. Off diagonal blocks are the computed
%	according to specified calling functions in the structure.
%	
%	The kernel structure was originally desined for operating with
%	single input motif kernels in gene networks, where correlations
%	between genes occur that are induced by an unobserved
%	transcription factor concentration.
%	
%	Information about the functions to call to compute the off
%	diagonal blocks is stored in a cell array called kern.block. By
%	default it is assumed that an off diagonal bloc will be computed
%	by a function called name1Xname2KernCompute. For example,
%	simXrbfKernCompute. If this function name doesn't exist
%	
%
%	KERN = MULTIKERNPARAMINIT(KERN) initialises the multiple output
%	block kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT, SIMXRBFKERNCOMPUTE


%	Copyright (c) 2006 Neil D. Lawrence



kern.nParams = 0;
kern.transforms = [];
if ~isfield(kern, 'comp')
  kern.comp=cell(0);
end

% the number of blocks should be the length of comp.
kern.numBlocks = length(kern.comp);
kern.isStationary = true;
for i = 1:length(kern.comp)
  if ~kern.comp{i}.isStationary
    kern.isStationary = false;
  end
  kern.comp{i} = kernParamInit(kern.comp{i});
  kern.nParams = kern.nParams + kern.comp{i}.nParams;
  kern.comp{i}.index = [];
  for j = 1:i-1
    func = [kern.comp{i}.type 'X' kern.comp{j}.type];
    if exist([func 'KernCompute']) == 2
      kern.block{i}.cross{j} = func;
      kern.block{i}.transpose(j) = false;
    else
      func = [kern.comp{j}.type 'X' kern.comp{i}.type];
      if exist([func 'KernCompute']) == 2
        kern.block{i}.cross{j} = func;
        kern.block{i}.transpose(j) = true;
      else
        warning('multiKernParamInit:noCrossKernel',['No cross covariance found between ' kern.comp{i}.type ...
                 ' and ' kern.comp{j}.type ' assuming independence.'])
        kern.block{i}.cross{j} = [];
        kern.block{i}.transpose(j) = 0;
      end
    end
  end
end
kern.paramGroups = speye(kern.nParams);

if isfield(kern, 'options') && isfield(kern.options, 'fixBlocks'),
  try,
    uuid = char(java.util.UUID.randomUUID);
    uuid = ['a' uuid(~(uuid=='-'))];
  catch,
    uuid = sprintf('%.30f', now);
    uuid = ['a' uuid(~(uuid=='.'))];
  end

  kern.fixedBlocks = zeros(1, kern.numBlocks);
  kern.fixedBlocks(options.fixBlocks) = 1;
  kern.uuid = uuid;

  multiKernCacheBlock(kern);
end
