function K = multiKernCacheBlock(kern, fhandle, arg, i, j, X, X2)

% MULTIKERNCACHEBLOCK
%
%	Description:
%
%	K = MULTIKERNCACHEBLOCK(KERN, FHANDLE, ARG, I, J, X1, X2) Maintains
%	a cache of values of multiKern blocks with constant parameters and
%	returns values
%	 Returns:
%	  K - the kernel matrix for the given inputs.
%	 Arguments:
%	  KERN - the structure containing the kernel.
%	  FHANDLE - function handle for computing the kernel
%	  ARG - cell array of arguments to fhandle
%	  I - the row of the block of the kernel to be computed.
%	  J - the column of the block of the kernel to be computed.
%	  X1 - first set of kernel inputs.
%	  X2 - second set of kernel inputs.
%
%	MULTIKERNCACHEBLOCK Clears the cache.
%	multiKernFixBlocks
%	
%
%	See also
%	MULTIKERNCREATE, MULTIKERNCOMPUTE, MULTIKERNCOMPUTEBLOCK, 


%	Copyright (c) 2008 Antti Honkela


persistent cache

if nargin == 1,
  cache.(kern.uuid) = cell(kern.numBlocks);
  K = [];
  return;
end

mycache = cache.(kern.uuid){i,j};

key = [X(:); X2(:)];

for k=1:length(mycache),
  if length(key) == length(mycache{k}{1}) && all(key == mycache{k}{1})
    K = mycache{k}{2};
    return;
  end
end

% No match if we get all the way here
if isempty(X2)
  K = fhandle(arg{:}, X);
else
  K = fhandle(arg{:}, X, X2);
end
cache.(kern.uuid){i, j}{end+1} = {key, K};
