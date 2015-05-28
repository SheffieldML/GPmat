function kern = multiKernFixBlocks(kern, blocks)

% MULTIKERNFIXBLOCKS
%
%	Description:
%
%	KERN = MULTIKERNFIXBLOCKS(KERN, BLOCKS) marks a number of blocks
%	fixed to allow caching their values
%	 Returns:
%	  KERN - the structure containing the kernel.
%	 Arguments:
%	  KERN - the structure containing the kernel.
%	  BLOCKS - vector of blocks to fix
%	
%
%	See also
%	MULTIKERNCREATE, MULTIKERNCACHEBLOCK


%	Copyright (c) 2008 Antti Honkela


try,
  uuid = char(java.util.UUID.randomUUID);
  uuid = ['a' uuid(~(uuid=='-'))];
catch,
  uuid = sprintf('%.30f', now);
  uuid = ['a' uuid(~(uuid=='.'))];
end

kern.fixedBlocks = zeros(1, kern.numBlocks);
kern.fixedBlocks(blocks) = 1;
kern.uuid = uuid;

multiKernCacheBlock(kern);
